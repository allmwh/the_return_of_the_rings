from pathlib import Path
import pandas as pd


def get_subset(human_df,subset):
    '''
    get subset(rbp, mrbp) from identified human_df,
    
    human_df: pd.DataFrame, 
    subset: 'rbp' 'mrbp'
    '''
    if subset == 'rbp':
        subset_list_path = Path("./rawdata/rbp_uniprotid_list.txt")
    elif subset == 'mrbp':
        subset_list_path = Path("./rawdata/mrbp_uniprotid_list.txt")
    else:
        raise ValueError("subset must be 'rbp' or 'mrbp'")
    
    with open(subset_list_path, "r") as tf:
        subset_list = tf.read().split('\n')
        
        
    subset_df = human_df[human_df['uniprot_id'].isin(subset_list)].reset_index(drop=True)
        
    return subset_df