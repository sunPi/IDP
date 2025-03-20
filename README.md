# IDP

Immune Deconvolution Pipeline - Deconvolutes TPM-normalized RNA-seq data by 
  methods such as 
  - 'quantiseq'
  - 'epic'
  - 'mcp_counter' 
  - 'abis' 
  - 'xcell'
  - 'consensus_tme' 
  - 'timer'

Usage: runIDP.R [options]

Options:
  -h --help                     Show this screen.
  -d --data_path=<PATH>         Path to the cleaned and deduped counts data set (RNA-seq data).
  -o --outfolder=<PATH>         Folder where the results are saved.
  -v --verbose=<BOOLEAN>        If set to 1, prints messages verbously.
  -V --version
