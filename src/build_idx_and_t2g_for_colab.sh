conda activate kbtools
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz

kallisto index gencode.v35.transcripts.fa.gz -i gencode.v35.transcripts.idx
kallisto index gencode.vM25.transcripts.fa.gz -i gencode.vM25.transcripts.idx


zgrep "^>" gencode.v35.transcripts.fa.gz | sed 's/>//g' | awk 'BEGIN {OFS = "\t"; FS = "|"}; {print $0, $2, $2}'  > v35.tr2g.tsv
# remove .\d+ ending
cat v35.tr2g.tsv |  awk -F'\t' -vOFS='\t' '{ gsub("\\.[0-9]*", "", $2) ; print }' | awk -F'\t' -vOFS='\t' '{ gsub("\\.[0-9]*", "", $3) ; print }' > v35.tr2gX.tsv

# mouse time
zgrep "^>" gencode.vM25.transcripts.fa.gz | sed 's/>//g' | awk 'BEGIN {OFS = "\t"; FS = "|"}; {print $0, $2, $2}'  > vM25.tr2g.tsv
# remove .\d+ ending
cat vM25.tr2g.tsv |  awk -F'\t' -vOFS='\t' '{ gsub("\\.[0-9]*", "", $2) ; print }' | awk -F'\t' -vOFS='\t' '{ gsub("\\.[0-9]*", "", $3) ; print }' > vM25.tr2gX.tsv

# Rscript bit to change the ENSG name from mouse to human 
# as the scVI scEiaD model uses the ENGS (human) gene names
module load R
Rscript ~/git/scEiaD/src/tr2g_humamize.R
# above outputs vM25.tr2gX.humanized.tsv

