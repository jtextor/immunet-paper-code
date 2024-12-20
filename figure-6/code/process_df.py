import pandas as pd
from scipy.spatial import KDTree
from config import slides_data_folder

activ_th = 0.4
mid_cd45ro_th = 0.33
high_cd45ro_th = 0.66

def get_neighbor_df(df):
    df_other = df.copy()

    df_other = df_other[df_other.celltype.isin(["CD3", "CD20", "CD4", "CD8", "FOXP3"])]

    df_other.loc[:,'cd45ro_expression'] = "low"
    df_other.loc[df_other.cd45ro > mid_cd45ro_th,'cd45ro_expression'] = "mid"
    df_other.loc[df_other.cd45ro > high_cd45ro_th,'cd45ro_expression'] = "hi"

    df_other = df_other.reset_index()
    tree = KDTree(df_other[['x','y']])
    neighbors = tree.query_pairs(r=14)
    for celltype in ["CD8", "CD4", "FOXP3", "CD20", "CD45RO"]:
        df_other.loc[:,celltype+"_neigh"] = 0

    for i,j in neighbors:
        type_i = df_other.loc[i, "celltype"]
        type_j = df_other.loc[j, "celltype"]

        df_other.loc[i, type_j+"_neigh"] += 1
        df_other.loc[j, type_i+"_neigh"] += 1

        df_other.loc[i, 'CD45RO_neigh'] += df_other.loc[j, 'cd45ro']
        df_other.loc[j, 'CD45RO_neigh'] += df_other.loc[i, 'cd45ro']
    return df_other


datasets = [f.name for f in slides_data_folder.iterdir() if f.is_dir()]
for dataset in datasets:
    dataset_folder = slides_data_folder / dataset
    slide_files = [i for i in dataset_folder.glob("*.csv") if "neighbors" not in str(i)]
    for file in slide_files:
        df = pd.read_csv(file)
        neighbor_df = get_neighbor_df(df)
        neighbor_df.to_csv(dataset_folder / (file.stem + ".neighbors.csv"), index=False)
