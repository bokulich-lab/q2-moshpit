# q2-moshpit
![CI](https://github.com/bokulich-lab/q2-moshpit/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-moshpit/branch/main/graph/badge.svg?token=PSCAYJUP01)](https://codecov.io/gh/bokulich-lab/q2-moshpit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**MO**dular **SH**otgun metagenome **P**ipelines with **I**ntegrated provenance **T**racking
QIIME 2 plugin for functional annotation and taxonomic classification of shotgun metagenomes.

## Installation
To install _q2-moshpit_, follow the installation steps described below.

```shell
mamba create -yn q2-shotgun \
  -c https://packages.qiime2.org/qiime2/2023.5/tested \
  -c bioconda -c conda-forge -c default q2-moshpit q2cli

conda activate q2-shotgun
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```

## Functionality
This QIIME 2 plugin contains actions used to annotate and classify (meta)genomes.
Below you will find an overview of actions available in the plugin.

| Action               | Description                                                | Underlying tool                                                    |
|----------------------|------------------------------------------------------------|--------------------------------------------------------------------|
| bin-contigs-metabat  | Bin contigs into MAGs using MetaBat 2.                     | [MetaBat 2](https://bitbucket.org/berkeleylab/metabat/src/master/) |
| build-kraken-db      | Fetch an existing or build a custom Kraken 2 database.     | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)                  |
| classify-kraken      | Classify reads/MAGs using Kraken 2.                        | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)                  |
| classify-kraken-bracken      | Classify reads using Kraken 2 and re-estimate abundances with Bracken.  | [Bracken](https://ccb.jhu.edu/software/bracken/index.shtml) |


## Dev environment
This repository follows the _black_ code style. To make the development slightly easier
there are a couple of pre-commit hooks included here that will ensure that your changes
follow that formatting style. Before you start working on the code, please
install the hooks by executing `make dev` in your conda environment. From then on,
they will be run automatically every time you commit any changes.