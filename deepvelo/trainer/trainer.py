from distutils.command.config import config
import numpy as np
from numpy import inf
import torch
import dgl

# from torchvision.utils import make_grid
from deepvelo.base import BaseTrainer
from deepvelo.utils import validate_config, inf_loop, MetricTracker
from deepvelo.logger import TensorboardWriter


class Trainer(BaseTrainer):
    """
    Trainer class
    """

    def __init__(
        self,
        model,
        criterion,
        metric_ftns,
        optimizer,
        config,
        data_loader,
        valid_data_loader=None,
        lr_scheduler=None,
        len_epoch=None,
    ):
        super().__init__(model, criterion, metric_ftns, optimizer, config)
        self.config = validate_config(config)
        self.data_loader = data_loader
        if len_epoch is None:
            # epoch-based training
            self.len_epoch = len(self.data_loader)
        else:
            # iteration-based training
            self.data_loader = inf_loop(data_loader)
            self.len_epoch = len_epoch
        self.valid_data_loader = valid_data_loader
        self.do_validation = self.valid_data_loader is not None
        self.lr_scheduler = lr_scheduler
        self.log_step = int(np.sqrt(data_loader.batch_size))

        self.train_metrics = MetricTracker(
            "loss", *[m.__name__ for m in self.metric_ftns], writer=self.writer
        )
        self.valid_metrics = MetricTracker(
            "loss", *[m.__name__ for m in self.metric_ftns], writer=self.writer
        )

    def _compute_core(self, batch_data):
        if isinstance(batch_data, dgl.nodeflow.NodeFlow):
            nf = batch_data
            nf.copy_from_parent()
            nf.layers[0].data["Ux_sz"] = nf.layers[0].data["Ux_sz"].to(self.device)
            nf.layers[0].data["Sx_sz"] = nf.layers[0].data["Sx_sz"].to(self.device)
            nf.layers[-1].data["Ux_sz"] = nf.layers[-1].data["Ux_sz"].to(self.device)
            nf.layers[-1].data["Sx_sz"] = nf.layers[-1].data["Sx_sz"].to(self.device)
            target = nf.layers[-1].data["velo"].to(self.device)
            output = self.model(nf)
        else:
            data_dict = batch_data
            x_u, x_s, target = data_dict["Ux_sz"], data_dict["Sx_sz"], data_dict["velo"]
            x_u, x_s, target = (
                x_u.to(self.device),
                x_s.to(self.device),
                target.to(self.device),
            )

            if self.config["arch"]["args"]["pred_unspliced"]:
                target_u = data_dict["velo_u"]
                target_u = target_u.to(self.device)
                # concate target to (batch, 2*genes), be careful of the order
                target = torch.cat([target, target_u], dim=1)

            output = self.model(x_u, x_s)
        return output, target

    def _smooth_constraint_step(self):
        topG = self.config["data_loader"]["args"]["topG"]
        neighbor_batch_ind = self.data_loader.dataset.gen_neighbor_batch(
            size=int(self.config["data_loader"]["args"]["batch_size"] / topG)
        )
        x_u = self.data_loader.dataset.Ux_sz[neighbor_batch_ind].to(
            self.device
        )  # (batch, genes)
        x_s = self.data_loader.dataset.Sx_sz[neighbor_batch_ind].to(self.device)
        output_neighbor = self.model(x_u, x_s)  # (batch, genes)
        self.optimizer.zero_grad()
        output_neighbor = output_neighbor.t().reshape([-1, topG])
        label_neighors = output_neighbor[:, 1:].detach()
        pivot = output_neighbor[:, 0:1]
        loss_c = torch.mean((label_neighors - pivot) ** 2)
        loss_c.backward()
        self.optimizer.step()

    def _train_epoch(self, epoch):
        """
        Training logic for an epoch

        :param epoch: Integer, current training epoch.
        :return: A log that contains average loss and metric in this epoch.
        """
        self.model.train()
        self.train_metrics.reset()
        if "mle" in self.config["loss"]["type"]:
            if "t+1" not in self.config["data_loader"]["args"]["type"]:
                if epoch <= self.config["trainer"]["guided_epochs"]:
                    self.data_loader.dataset.neighbor_time = 1
                    print(
                        f"This is a guided_epoch, use neighbors at t{self.data_loader.dataset.neighbor_time}."
                    )
                else:
                    self.data_loader.dataset.neighbor_time = 0
        if (not self.data_loader.shuffle) and self.data_loader.is_large_batch:
            loader = self.data_loader.dataset.large_batch(self.device)
        else:
            loader = self.data_loader
        for batch_idx, batch_data in enumerate(loader):
            output, target = self._compute_core(batch_data)
            if self.config["mask_zeros"]:
                data_dict = batch_data
                mask = data_dict["mask"].to(self.device)  # (batch, n_gene)
                output = output * mask
                target = target * mask
            self.optimizer.zero_grad()
            loss = (
                self.criterion(output, target)
                if "mle" not in self.config["loss"]["type"]
                else self.criterion(
                    output,
                    current_state=(
                        torch.cat(
                            [
                                batch_data["Sx_sz"],
                                batch_data["Ux_sz"],
                            ],
                            dim=1,
                        ).to(self.device)
                        if self.config["arch"]["args"]["pred_unspliced"]
                        else batch_data["Sx_sz"].to(self.device)
                    ),
                    idx=batch_data["t+1 neighbor idx"].to(self.device),
                    candidate_states=self.candidate_states,
                    spliced_counts=batch_data["Sx_sz"].to(self.device),
                    unspliced_counts=batch_data["Ux_sz"].to(self.device),
                    **self.config["loss"]["args"],
                )
            )
            loss.backward()
            self.optimizer.step()

            if self.config["constraint_loss"]:
                self._smooth_constraint_step()

            self.writer.set_step((epoch - 1) * self.len_epoch + batch_idx)
            self.train_metrics.update("loss", loss.item())
            for met in self.metric_ftns:
                self.train_metrics.update(met.__name__, met(output, target))

            if batch_idx % self.log_step == 0:
                self.logger.debug(
                    "Train Epoch: {} {} Loss: {:.6f}".format(
                        epoch, self._progress(batch_idx), loss.item()
                    )
                )
                # self.writer.add_image('input', make_grid(data.cpu(), nrow=8, normalize=True))

            if batch_idx == self.len_epoch:
                break
        log = self.train_metrics.result()

        if self.do_validation:
            val_log = self._valid_epoch(epoch)
            log.update(**{"val_" + k: v for k, v in val_log.items()})

        if self.lr_scheduler is not None:
            self.lr_scheduler.step()
        return log

    def _valid_epoch(self, epoch):
        """
        Validate after training an epoch

        :param epoch: Integer, current training epoch.
        :return: A log that contains information about validation
        """
        self.model.eval()
        self.valid_metrics.reset()
        with torch.no_grad():
            for batch_idx, batch_data in enumerate(self.valid_data_loader):
                output, target = self._compute_core(batch_data)
                loss = self.criterion(output, target)

                self.writer.set_step(
                    (epoch - 1) * len(self.valid_data_loader) + batch_idx, "valid"
                )
                self.valid_metrics.update("loss", loss.item())
                for met in self.metric_ftns:
                    self.valid_metrics.update(met.__name__, met(output, target))
                # self.writer.add_image('input', make_grid(data.cpu(), nrow=8, normalize=True))

        # add histogram of model parameters to the tensorboard
        for name, p in self.model.named_parameters():
            self.writer.add_histogram(name, p, bins="auto")
        return self.valid_metrics.result()

    def _progress(self, batch_idx):
        base = "[{}/{} ({:.0f}%)]"
        if hasattr(self.data_loader, "n_samples"):
            current = batch_idx * self.data_loader.batch_size
            total = self.data_loader.n_samples
        else:
            current = batch_idx
            total = self.len_epoch
        return base.format(current, total, 100.0 * current / total)

    def eval(self, eval_loader, return_kinetic_rates=False):
        """
        Evaluate the model on a given dataset, provided by eval_loader

        :param model: model to evaluate
        :param eval_loader: dataset loader containing dataset to evaluate on
        """
        self.model.eval()
        n_genes = self.config["arch"]["args"]["n_genes"]
        velo_mat = []
        velo_mat_u = []
        kinetic_rates = {}
        if (not eval_loader.shuffle) and eval_loader.is_large_batch:
            loader = eval_loader.dataset.large_batch(self.device)
        else:
            loader = eval_loader
        with torch.no_grad():
            for batch_idx, batch_data in enumerate(loader):
                output, _ = self._compute_core(batch_data)
                pred_s = output[:, 0:n_genes]
                velo_mat.append(pred_s.cpu().data)
                if self.config["arch"]["args"]["pred_unspliced"]:
                    pred_u = output[:, n_genes:]
                    velo_mat_u.append(pred_u.cpu().data)

                # record kinetic rates
                if return_kinetic_rates:
                    cur_kinetic_rates = self.model.get_current_batch_kinetic_rates()
                    for k, v in cur_kinetic_rates.items():
                        if k not in kinetic_rates:
                            kinetic_rates[k] = []
                        kinetic_rates[k].append(v.cpu().data)

            velo_mat = np.concatenate(velo_mat, axis=0)
            if return_kinetic_rates:
                for k, v in kinetic_rates.items():
                    kinetic_rates[k] = np.concatenate(v, axis=0)
            if self.config["arch"]["args"]["pred_unspliced"]:
                velo_mat_u = np.concatenate(velo_mat_u, axis=0)
        return velo_mat, velo_mat_u, kinetic_rates
