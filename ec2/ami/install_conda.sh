wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
