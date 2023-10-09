
sudo apt update
sudo apt install -y software-properties-common

sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer

# install conda in /opt/miniconda3

# install mamba as root:
sudo su - 
conda install -n base -c conda-forge mamba

# install snake_course environment
mamba env create -f conda/snakemake_course.yaml

# install snakemake
conda activate snake_course
mamba install -c conda-forge -c bioconda snakemake
