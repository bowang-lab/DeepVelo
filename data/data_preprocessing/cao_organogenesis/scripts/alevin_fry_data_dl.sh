#!/bin/bash 

mkdir -p data/cao_organogenesis
cd data/cao_organogenesis

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
wget https://zenodo.org/record/5676291/files/cell_annotations.txt
wget https://zenodo.org/record/5676291/files/transcriptome_splici_fl52_t2g_3col.tsv
wget https://zenodo.org/record/5676291/files/cell_barcodes.txt
wget https://zenodo.org/record/5676291/files/gene_annotation_table.txt
wget https://zenodo.org/record/5676291/files/splici_idx_mouse_gencodeM25.tar.gz
