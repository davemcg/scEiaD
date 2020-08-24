import sys
import parc
import pandas as pd

res = float(sys.argv[2])
reduction = pd.read_csv(sys.argv[1], index_col  = 0)
parc1 = parc.PARC(reduction, jac_std_global= 'median', random_seed =1, small_pop = 100, resolution_parameter = res)
parc1.run_PARC()

parc_labels = parc1.labels
out = pd.DataFrame( pd.Categorical(parc_labels))
out['Barcode'] = reduction.index
out.to_csv(sys.argv[3])

