import re
import pandas as pd
from Bio import SeqIO
from pathlib import Path

def get_protein_name(uniprot_id,human_df):
    '''
    輸入uniprot_id拿到protein的名字
    '''
    protein_name = human_df[human_df['uniprot_id']==uniprot_id]['protein_name'].values[0]
    gene_name = human_df[human_df['uniprot_id']==uniprot_id]['gene_name'].values[0]
    
    #判斷有括號
    if "(" in protein_name:
        protein_name = re.findall( r"[\w\s.,'-]*", string = protein_name )
        protein_name = protein_name[0]
    else:
        protein_name = protein_name
    
    protein_name = re.sub("/","_",protein_name)
    protein_name = re.sub(" ","_",protein_name)
    
    return {"protein_name":protein_name,
            "gene_name":gene_name}

def get_uniprot_rawdata(path):
    '''
    path: str, the uniprot tab file path
    
    Load protein data downloaded from uniprot (https://www.uniprot.org/uploadlists/), 
    the columns of uniprot must be 
    "Entry, Gene names (primary), Protein names, Sequence, Organism ID" 
    and save as "Tab-separated format (.tab)" 
    '''
    df = pd.read_csv(path,sep='\t',names=['uniprot_id','gene_name','protein_name','protein_sequence','taxonomy'])
    df = df.drop(0).reset_index().drop(axis=1,labels='index')
    return df

def get_uniprot_id_from_fasta(path):
    '''
    從fastapath拿到uniprotid
    '''
    path = Path(path)
    uniprot_id = path.parts[-1].split('.')[0]
    return uniprot_id

def fasta_to_seqlist(path):
    '''
    輸入fasta讀成list

    path: fasta檔案

    return: list
    '''
    return list(SeqIO.parse(str(path),'fasta'))

def find_human_sequence(path):
    '''
    輸入fasta檔案，找到人類那條

    path: 要找的fasta檔案
    
    return: 人類的序列
    '''
    
    path = Path(path)

    #讀檔
    all_sequence_list = fasta_to_seqlist(path)
    uniprot_id = get_uniprot_id_from_fasta(path)

    #找出human那條
    for i in all_sequence_list:
        tax_id = i.description.split("|")[2]
        if int(tax_id) == 9606:
            sequence = str(i.seq)
            return {"uniprot_id":uniprot_id,
                    "sequence":sequence}
    
    #沒human的err handle
    print("{}, FASTA PATH {} DOES NOT HAVE HUMAN SEQUENCE".format(uniprot_id,str(path)))
    raise Exception("error")




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
        raise ValueError("subset must be 'rbp' or 'mrbp' ")
    
    with open(subset_list_path, "r") as tf:
        subset_list = tf.read().split('\n')
        
        
    subset_df = human_df[human_df['uniprot_id'].isin(subset_list)].reset_index(drop=True)
        
    return subset_df


def get_taxid_dict():
        '''
        物種id，之後可以創更多
        '''
        return {
                131567:"cellular_organisms",
                2759:"eukaryota",
                33154:"opisthokonta",
                33208:"metazoa",
                6072:"eumetazoa",
                33213:"bilateria",
                33511:"deuterostomia",
                7711:"chordata",
                89593:"craniata",
                7742:"vertebrata",
                7776:"gnathostomata",
                117570:"teleostomi",
                117571:"euteleostomi",
                8287:"sarcopterygii",
                1338369:"dipnotetrapodomorpha",
                32523:"tetrapoda",
                32524:"amniota",
                40674:"mammalia",
                32525:"theria",
                9347:"eutheria",
                1437010:"boreoeutheria",
                314146:"euarchontoglires",
                9443:"primates",
                376913:"haplorrhini",
                314293:"simiiformes",
                9526:"catarrhini",
                314295:"hominoidea",
                9604:"hominidae",
                207598:"homininae",
                9605:"homo"
               }