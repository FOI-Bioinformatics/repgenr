# RepGenR

## Introduction
RepGenR (Representative-Genome Repositories) is a bioinformatics tool designed to streamline the process of clustering genome sequences. It uses average nucleotide identity (ANI) to produce representative genome repositories, providing researchers with a more manageable dataset for large-scale genomic studies.

## Features
RepGenR offers a suite of features to facilitate the organization and analysis of genomic data:

- Taxa selection from the GTDB (Genome Taxonomy DataBase).
- Downloading and formatting genome sequences into a standardized naming convention.
- De-replication of genomes based on average nucleotide identity (ANI).
- Phylogenetic tree computation for both dereplicated and all sequence datasets.
- Generation of phylogenetic trees in a format compatible with downstream processes such as FlexTaxD.

## Getting Started
To begin using RepGenR, please consult the [wiki](https://github.com/FOI-Bioinformatics/repgenr/wiki) for detailed installation instructions and usage guidelines.

## Installation
Quick installation guide for users who are familiar with the command line:

```bash

# Create a new environment and install dependencies. Make sure to activate the environment before installing RepGenR.
mamba create -n RepGenR -y python=3 matplotlib drep checkm-genome ncbi-datasets-cli mashtree progressivemauve iqtree ete3

# Clone the repository
git clone https://github.com/FOI-Bioinformatics/repgenr.git

# Navigate to the RepGenR directory
cd repgenr
pip install .

# Follow the detailed installation instructions in the wiki
```
## Usage
For a complete guide on how to use RepGenR, please refer to the [wiki](https://github.com/FOI-Bioinformatics/repgenr/wiki). Here's a quick start command to get you running:

```bash
# Replace with actual usage command
repgenr --help
```
See [module wiki page](https://github.com/FOI-Bioinformatics/repgenr/wiki/2.-Modules) for available sub commands.

## Contributing
We welcome contributions to RepGenR! If you have suggestions or contributions, please open an issue or pull request.

## License
RepGenR is licensed under the [MIT License](LICENSE).
