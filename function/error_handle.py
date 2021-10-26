from pathlib import Path

def no_human_sequence(uniprot_id,path):
    '''
    for utilities.find_human_sequence use
    看find_human_sequence裡面沒找到9606丟錯誤
    '''

    print("{}, FASTA PATH {} DOES NOT HAVE HUMAN SEQUENCE".format(uniprot_id,str(path)))
    raise Exception("error")

def no_fasta_file(uniprot_id,path):
    '''
    沒有fasta檔案
    '''
    path = Path(path)
    if not path.is_file():
            print("{} DOES NOT HAVE FASTA FILE".format(uniprot_id))
            raise Exception("error")


def od_ident_score_diff_length(uniprot_id,score,od_ident):
    '''
    od_ident:是uniprot資料庫的長度
    score:算entropy的長度

    因為oma資料庫的關係，長度會不一樣
    '''
    if len(score) != len(od_ident):
        print("{} ENTROPY LENGTH IS NOT EQUAL WITH OD_IODENT".format(uniprot_id))
        raise Exception("error")



