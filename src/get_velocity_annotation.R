library(tidyverse)
library(BUSpaRse)
args = commandArgs(trailingOnly = T)
gtf <- args[1]
ref <- args[2]
tech <- args[3]
outdir <- args[4]
    
n=0
tech2rl <- case_when(tech == '10xv2' ~ 98, 
                     tech == '10xv3' ~ 91, 
                     tech == "SMARTSeq_v2" ~ 50, 
                     tech == "DropSeq"  ~ 50, # not super sure about this one but fux it 
                     tech == "SMARTerSeq_v3" ~ 50, 
                     tech == "SCRBSeq"  ~ 100, 
                     tech == "SMARTSeq_v4" ~ 75,
                     tech == "C1" ~ 94)


if(ref == 'hs-homo_sapiens'| ref == 'DNTX'){
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

} else{
    message('ERROR')
    quit()
}