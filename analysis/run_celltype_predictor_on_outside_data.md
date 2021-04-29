# Setup new conda environment 
```
conda create -n scEiaD_CT_learner
conda activate scEiaD_CT_learner
conda install python=3.8 pandas numpy scikit-learn xgboost=1.3
pip install --quiet kb-python
pip install --quiet scvi-tools[tutorials]==0.9.0
```
# download scVI model which will project your data into the the batch corrected scEiaD space
```
# ~13mb
wget -O scVI_scEiaD.tgz https://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/2021_03_17__scVI_scEiaD.tgz
tar -xzf scVI_scEiaD.tgz
```
# download xgboost data
```
## ~12 mb
wget -O celltype_ML_model.tar https://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_03_17/2021_cell_type_ML_all.tar
tar -xzf celltype_ML_model.tar
```
# download two python scripts
```
## the first contain the xgboost ML functions
wget -O celltype_predictor.py https://raw.githubusercontent.com/davemcg/scEiaD/master/src/cell_type_predictor.py
## the second is a simple workflow that runs scVI to get the latent dims then runs the xgboost ML prediction
wget -O scEiaD_celltype_predictor.py https://raw.githubusercontent.com/davemcg/scEiaD/master/src/scEiaD_celltype_predictor.py
```

# get the human (or mouse) kallisto index and transcript to gene tsv 
```
# ~3 GB
# human
wget https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/gencode.v35.transcripts.idx
wget https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/v35.tr2gX.tsv
# mouse
wget https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/gencode.vM25.transcripts.idx
wget https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/vM25.tr2gX.tsv
```
# run kb count
```
## (tweak the tech as needed - you may not be using 10xv2)
## kb count --overwrite --h5ad -i gencode.v35.transcripts.idx -g v35.tr2gX.tsv -x DropSeq -o output_dir --filter bustools -t 12 \
## sample_R1.fastq.gz sample_R2.fastq.gz
# toy example with human retinal organoids
wget -O sample_1.fastq.gz https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/SRR11799731_1.head.fastq.gz
wget -O sample_2.fastq.gz https://hpc.nih.gov/~mcgaugheyd/scEiaD/colab/SRR11799731_2.head.fastq.gz

kb count --overwrite --h5ad -i gencode.v35.transcripts.idx -g v35.tr2gX.tsv -x DropSeq -o output_dir --filter bustools -t 4 sample_1.fastq.gz sample_2.fastq.gz
```
# do cell type prediction!
```
## args:
## 1. adata.h5ad 
## 2. output csv which will contain the probability that each cell (row) is a cell type (column) 
## the last column is the most likely cell type. 
## 3. number of epochs that scVI should run for (more epochs is more time to tune the model error rate)
## The scVI people recommend around 50. We have found that far fewer (5!) also works quite well
## 4. cutoff for cell type calling failure. if the highest prob is <USER SPECIFIED THRESHOLD (e.g. 0.5 is 50%) 
## then the cell type will be called "None"
## 5. "mouse" gets special handling to deal with gene name issues (as scEiaD uses human ensembl gene IDs)
python scEiaD_celltype_predictor.py output_dir/counts_filtered/adata.h5ad celltypes.csv 50 0.5 mouse
# get number of columns and then get quick stats
awk -F"," '{print NF}' celltypes.csv | uniq
# the answer is 32
cut -d, -f32 celltypes.csv | sort | uniq -c
```
