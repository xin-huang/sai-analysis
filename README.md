[![license](https://img.shields.io/badge/license-GPL%20v3-black.svg?style=flat-square)](LICENSE)
[![build Status](https://img.shields.io/github/actions/workflow/status/xin-huang/sai-analysis/dry-run.yaml?branch=main&style=flat-square&label=dry-run)](https://github.com/xin-huang/sai-analysis/actions)

# sai-analysis

## Introduction

This repository contains a Snakemake workflow designed to reproduce the results of a scan for candidates of adaptive introgression using [sai](https://github.com/xin-huang/sai). The workflow has been tested on Oracle Linux 9 using the Life Science Compute Cluster at the University of Vienna.

## Usage

1. Install [Mambaforge (version: 23.3.1)](https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Mambaforge-23.3.1-1-Linux-x86_64.sh).

2. Clone this repository:

```
git clone https://github.com/xin-huang/sai-analysis
cd sai-analysis
```

3. Create the environment:

```
mamba env create -f workflow/envs/env.yaml
```

4. Activate the environment:

```
mamba activate sai-analysis
```

5. Run the analysis locally:

```
snakemake -c 1
```

6. Run the analysis on HPC:

```
snakemake -c 1 --profile config/slurm
```

Users should adjust the resource parameters in each Snakemake file to match their cluster settings and modify the `config.yaml` file in `config/slurm` to suit their job scheduler.
