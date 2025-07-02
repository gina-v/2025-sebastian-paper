# Workflow to help Sebastian with Adherent  E. Coli and S. Typhiomerium

Here we will annotate the workflow steps that must be done to create a group of contigs for blasting purposes
 
But, before we do anything lets look at conda....

go to your `.condarc` file that stores your personal configuration for your conda install:

```
nano ~/.condarc
```

```
channels:
  - conda-forge
  - bioconda
  - defaults
envs_dirs:
  - /share/apps/conda/environments
pkgs_dirs:
  - /share/apps/conda/pkgs
  - ~/.conda/pkgs

channel_priority: strict
```

I started by downloading `megahit`, `fastqc`, and `fastp` into the `sebastian.yml` env.

The code below creates a yaml config file for the env name sebastian:
```
`conda env export --from-history --name sebastian > "sebastian_env.yml"`
```

To install:
```
conda env create --file sebastian_env.yml
```
