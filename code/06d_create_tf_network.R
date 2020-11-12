

# try using PANDA to create a TF network for the SAE data
library('pandaR')

data(pandaToyData)
# motif, ppi, expression
# I think we probably want "data/tissue_motif.txt":
#   this is 19mil edges big though