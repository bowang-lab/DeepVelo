# Analysis scripts for DeepVelo

This folder contains the code necessary to reproduce the analysis for the DeepVelo manuscript

### Instructions:

Most notebooks can run directly and download necessary data automatically. For running python notebooks using the hindbrain data or running the R scripts:

1. Please download the necessary files:

```
wget -O deepvelo_data.tar.gz https://www.dropbox.com/s/xap6r5bh6er93uh/deepvelo_data.tar.gz?dl=0
```

2. Untar and extract the main directory

```
tar -xzvf deepvelo_data.tar.gz
```

**For R scripts** - please change directories into the `r_analysis` subfolder and follow the instructions
to download and install the conda environment for generating the R figures and analyses.
