#!/bin/bash 
#scVI
conda create -y -c bioconda -c conda-forge -n scVI scvi
#scanorama
conda create -y -c conda-forge -n scanorama python=3.7 python-annoy 
conda activate scanorama
pip install 
conda deactivate
#magic 
conda create -y -n magic python=3.7
conda activate magic 
pip install magic-impute
conda deactivate
#bbknn 
conda create -y -c bioconda -n bbknn bbknn
#INSCT
conda create -y -n INSCT python
conda activate INSCT
pip install git+https://github.com/lkmklsmn/insct.git
conda deactivate
#DESC
conda create -y -n DESC python=3.6
conda activate DESC
pip install desc
conda deactivate
#parc 
conda create -y -n parc python=3.6 llvmlite pybind11 numpy 
conda activate parc 
pip install parc
conda deactivate
#phate
conda create -y -n phate python=3.6
conda activate phate
pip install phate 
conda deactivate
#scIB
conda create -y -n scIB python=3.6
conda activate scIB
pip install git+https://github.com/theislab/scib.git
conda deactivate
#integrate_scRNA
conda create -y -n integrate_scRNA python=3.6 pandas numpy xgboost scikit-learn
