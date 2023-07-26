# Single cell ATAC-seq data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155178
import pandas as pd

scatac_file = pd.read_csv('../GSE155178_ACR_x_cell.binary.sparse.txt.gz',sep='\t',header = None)

scatac_file.columns = ['loc','tag','num']
# clean the raw input data
df = scatac_file
s = pd.DataFrame(df['loc'].str.split("_").to_list())
df = df.join(s)
df.columns = ['loc','tag','num','chr','start','end']
df = df[['chr','start','end','tag','num']]

# get rid of "CB:Z:" and save the real cell barcode info
s = pd.DataFrame(df['tag'].str.split(":").to_list())
df = df.join(s)
df.columns = ['chr','s','e','frg','num','CB','Z','real_frg']
df = df[['chr','s','e','tag','num']]
# get the tissue info
s = pd.DataFrame(df['tag'].str.split("-").to_list())
df = df.join(s)
df.columns = ['chr','s','e','tag','num','real_tag','tissue']
df = df[['chr','s','e','real_tag','tissue','num']]
df['real_tag'] = df['real_tag'].str.replace('CB:Z:', '', regex=True)
# get rid of "chr" to match the B73 maize genome annotation.
df["chr"] = df["chr"].str.replace("chr", "")
# get the unique tissue list
# save single cell ATAC-seq data based on the tissue type.
tissue_list = list(df.tissue.unique())

for t in tissue_list:
    tissue_df = df[df['tissue']==t]
    tissue_df = tissue_df[['chr','s','e','real_tag','num']]
    print(t)
    print(tissue_df.head())
    tissue_df.to_csv('./scATAC/B73_tissue_Frg/GSE155178_'+t+'_real_frg.tsv.gz', sep="\t",index = False, header = False,compression='gzip')

print("DONE!!")
