import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

base_dir = "/home/mrama/Desktop/MD/eMASS-MD_complete_data/enzyme_models/ENO/ENO_param_scan2/output/treated_data/"

file_in = ''.join([base_dir, 'rateconst_ENO_customRatio_1_0.1_1.csv'])

print(file_in)

data_df = pd.read_csv(file_in, sep='\t')

sns.clustermap(data=data_df,  metric='euclidean', 
            col_cluster=False, center=0,
            vmin=-6, vmax=9, yticklabels=False, cmap='RdBu_r')

plt.savefig(file_in[:-4] + '.pdf')
plt.close()
