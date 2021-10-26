from tqdm.notebook import trange
import requests
import json
import re
from pathlib import Path 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from function.utilities import fasta_to_seqlist
from function.utilities import get_taxid_dict
from function.utilities import find_human_sequence

class FetchOmaSeq():
    def __init__(self):
        pass

    def get_oma_seq(self,uniprot_id,path):
        try:
            path = Path(path)
            fasta_path = path / "{}.fasta".format(uniprot_id)
            orthologs_list = self.__get_orthologs(uniprot_id)
            self.__get_fasta(orthologs_list,fasta_path)
        except:
            print("{} OMA DOES NOT HAVE THIS UNIPROT_ID RECORD".format(uniprot_id))
            raise Exception("error")

        #uniprot_id_consist_check，uid跟oma的不合check
        uniprot_id_oma_fassta = find_human_sequence(fasta_path)['uniprot_id']
        if uniprot_id != uniprot_id_oma_fassta:
            fasta_path.unlink()
            print("{} IS NOT CONSIST WITH OMA's record".format(uniprot_id))
            raise Exception("error")
                
    def pipeline_get_oma_seq(self,uniprot_ids,path):
        path = Path(path)
        t = trange(len(uniprot_ids), leave=True,position=0)
        for i in t:
            self.get_oma_seq(uniprot_ids[i],path)

    def __get_protein_info_from_entry(self,uniprot_id):
        resp = requests.get('https://omabrowser.org/api/protein/{}/'.format(uniprot_id))
        oma_raw = json.loads(resp.text)
        
        species = oma_raw['species']['species']
        #去掉可悲的括號
        species = re.sub("\(.*\)","",species) 

        oma_id = oma_raw['omaid']
        canonical_id = oma_raw['canonicalid']
        taxon_id = oma_raw['species']['taxon_id']
        sequence = oma_raw['sequence']
        
        return {'species':species,
                'oma_id':oma_id,
                'canonical_id':canonical_id,
                'taxon_id':taxon_id,
                'sequence':sequence}


    def __get_orthologs(self,uniprot_id):
        resp = requests.get('https://omabrowser.org/api/group/{}/'.format(uniprot_id))
        oma_raw = json.loads(resp.text)
        
        orthologs_list = []
        t = trange(len(oma_raw['members']), desc=uniprot_id, leave=False,position=2)
        
        for i in t:
            orthologs_list.append(self.__get_protein_info_from_entry(oma_raw['members'][i]['entry_nr'])) #拿data
        return orthologs_list 

    def __get_fasta(self,orthologs_list,path):
        fasta_list = []
        for i in orthologs_list:
            record = SeqRecord(Seq(i['sequence']),
                    id = i['oma_id'], 
                    description= "| {} | {} | {}".format(i['species'],i['taxon_id'],i['canonical_id']))
            fasta_list.append(record)
        SeqIO.write(fasta_list,path, "fasta")
    


class TaxSeqFilter():
    def __init__(self,taxonomy):
        '''
        輸入taxid，去過濾掉oma來的序列

        '''
        #species filter
        resp = requests.get('https://omabrowser.org/api/taxonomy/{}'.format(taxonomy))
        self.taxonomy = taxonomy
        self.taxonomy_list = resp.text
        
    def taxfilter(self,oma_fasta_path,grouped_fasta_path):
        '''
        輸入fasta_seq和輸出路徑，過濾fasta

        oma_fasta_path: Path,
        grouped_fasta_path: Path
        '''

        #read
        oma_fasta_list = fasta_to_seqlist(oma_fasta_path)
        
        #filter
        filtered_list = []
        for i in oma_fasta_list:
            tax_id = i.description.split("|")[2].replace(' ', '')
            
            if tax_id in self.taxonomy_list:
                filtered_list.append(i)
        
        with open(grouped_fasta_path, "w") as output_handle:
            SeqIO.write(filtered_list, output_handle, "fasta")

    def taxfilter_pipeline(self,oma_path,grouped_tax_path):
        '''
        直接給oma根目錄，就開始篩物種，放到要的資料夾

        oma_path: Path()
        grouped_tax_path: 目標資料夾
        '''
        
        tax_id = Path(str(self.taxonomy))
        species = get_taxid_dict()[self.taxonomy]
        grouped_tax_path.mkdir(exist_ok=True,parents=True)
        
        fasta_pathlist = list(Path(oma_path).rglob('*.fasta'))
        t = trange(len(fasta_pathlist), desc="", leave=True,position=2)

        for i in t:
            uniprot_id_path = fasta_pathlist[i]
            
            t.set_description("{} ({}): {}".format(species,tax_id,uniprot_id_path.parts[-1].split(".")[0]))
            self.taxfilter(uniprot_id_path,grouped_tax_path/Path(uniprot_id_path.parts[-1]))
    
