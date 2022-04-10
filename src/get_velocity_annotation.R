library(tidyverse)
library(BUSpaRse)
library(glue)
args = commandArgs(trailingOnly = T)
gtf <- args[1]
ref <- args[2]
tech <- args[3]
outdir <- args[4]
git_dir = args[5]
    
n=0
tech2rl <- case_when(tech == '10xv2' ~ 98, 
                     tech == '10xv3' ~ 91,  
                     tech == "DropSeq"  ~ 50, # not super sure about this one but fux it 
                     tech == 'well' ~ 100)



if (ref == 'hs-homo_sapiens'| ref == 'DNTX'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    get_velocity_files(gtf, L = tech2rl, Genome = BSgenome.Hsapiens.UCSC.hg38, 
                       out_path =outdir , 
                       isoform_action = "separate")
        
    
} else if (ref == 'mm-mus_musculus'){
    library(BSgenome.Mmusculus.UCSC.mm10)
    get_velocity_files(gtf, L = tech2rl, Genome = BSgenome.Mmusculus.UCSC.mm10, 
                       out_path =outdir , 
                       isoform_action = "separate")
    

} else if (ref == 'mf-macaca_mulatta') {
    library(BSgenome.Mmulatta.UCSC.rheMac10)
    get_velocity_files(gtf, L = tech2rl, Genome = BSgenome.Mmulatta.UCSC.rheMac10, 
                   out_path =outdir , 
                   isoform_action = "separate")

} else if (ref =='gg-gallus_gallus') {
	library(BSgenome.Ggallus.UCSC.galGal6)
    get_velocity_files(gtf, L = tech2rl, Genome = BSgenome.Ggallus.UCSC.galGal6, 
                   out_path =outdir , 
                   isoform_action = "separate")
} else{
    message('ERROR')
    quit()
}

## this is definetly not the best way to do this but it was the easiest with out changing a ton of code 
if(tech == 'well'){
    print(glue('python3 {git_dir}/src/remove_entry_from_fasta.py --infasta {outdir}/cDNA_introns.fa --txToRemove {outdir}/introns_tx_to_capture.txt  --outfasta {outdir}/cDNA_only.fa && rm {outdir}/cDNA_introns.fa && mv {outdir}/cDNA_only.fa  {outdir}/cDNA_introns.fa'))
    system(glue('python3 {git_dir}/src/remove_entry_from_fasta.py --infasta {outdir}/cDNA_introns.fa --txToRemove {outdir}/introns_tx_to_capture.txt  --outfasta {outdir}/cDNA_only.fa && rm {outdir}/cDNA_introns.fa && mv {outdir}/cDNA_only.fa  {outdir}/cDNA_introns.fa'))

}

