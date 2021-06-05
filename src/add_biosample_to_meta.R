# add biosample to old metadata
orig_meta <- read_tsv('~/git/scEiaD/data/sample_run_layout_organism_tech.tsv')
grab <- efetch('SRR13075732', db = 'sra', retmode = 'xml')
  str_extract(grab$content, 'SAMN\\d+')