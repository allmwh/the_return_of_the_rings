import re
import json
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from tqdm.notebook import trange
from Bio.SeqRecord import SeqRecord

from function.utilities import get_taxid_dict
from function.utilities import fasta_to_seqlist
from function.utilities import find_human_sequence


def __uniprot_id_consistance_check(fasta_path,uniprot_id):
    # some uniprot id in OMA paralogs is not consist with uniprot 
        uniprot_id_oma_fassta = find_human_sequence(fasta_path)["uniprot_id"]
        if uniprot_id != uniprot_id_oma_fassta:
            fasta_path.unlink()
            raise Exception("{} in uniprot is not consist with OMA's record, delete this record".format(uniprot_id))
class FetchOmaSeqBatch():
    '''
    抓fasta速度比較快的方法:
    從  https://omabrowser.org/oma/omagroup/Q13148/fasta/ 拿到他們格式的fasta，
    再從 https://omabrowser.org/api/group/Q13148/ 拿到fasta資料，然後改我自己的格式
    '''
    def __init__(self):
        pass

    def get_oma_seq(self, uniprot_id, path):
        '''
        pipeline: 從oma抓fasta, 拿到fasta詳細資料, 改fasta資訊跟原本我的一樣
        '''

        oma_path = Path(path)
        oma_fasta_path = oma_path / "{}.fasta".format(uniprot_id)

        # get raw fasta
        self.__get_oma_fasta(uniprot_id, oma_fasta_path)

        # get fasta info
        fasta_info_dict = self.__get_fasta_info(uniprot_id)

        # get mod info fasta
        self.__mod_fasta_info(oma_fasta_path, oma_fasta_path, fasta_info_dict)

        # uniprot id consistance check
        __uniprot_id_consistance_check(oma_fasta_path, uniprot_id)

    def __get_oma_fasta(self, uniprot_id, fasta_path):
        '''
        從oma拿到整包的fasta
        '''
        try:
            url = "https://omabrowser.org/oma/omagroup/{}/fasta/".format(uniprot_id)
            resp = requests.get(url)
            resp.raise_for_status()
            with open(fasta_path, "w") as file:
                file.write(resp.text)
        except:
            raise Exception("{} OMA fetch fasta failed".format(uniprot_id))

    def __get_fasta_info(self, uniprot_id):
        '''
        拿到所要的fasta資料
        '''
        try:
            url = "https://omabrowser.org/api/group/{}/".format(uniprot_id)
            resp = requests.get(url)
            resp.raise_for_status()
            oma_raw = json.loads(resp.text)
            fasta_info_dict = {}
            for i in oma_raw['members']:

                species = i["species"]["species"]
                # species name is too long, remove some strain info
                species = re.sub("\(.*\)", "", species)

                oma_id = i["omaid"]
                canonical_id = i["canonicalid"]
                taxon_id = i["species"]["taxon_id"]

                fasta_info_dict[oma_id] = {
                    "oma_id": oma_id,
                    "species": species,
                    "canonical_id": canonical_id,
                    "taxon_id": taxon_id,
                }

            return fasta_info_dict
        except:
            raise Exception("{} OMA fetch fasta seqeuence info failed".format(uniprot_id))

    def __mod_fasta_info(self, oma_fasta_path, mod_fasta_path, fasta_info_dict):
        '''
        用拿到的fasta資料去修改成新的fasta
        '''

        fasta_list = list(SeqIO.parse(str(oma_fasta_path), 'fasta'))
        mod_fasta_list = []
        for seq_record in fasta_list:
            id = seq_record.id
            record = SeqRecord(seq=seq_record.seq,
                               id=id,
                               description="| {} | {} | {}".format(fasta_info_dict[id]["species"],
                                                                   fasta_info_dict[id]["taxon_id"],
                                                                   fasta_info_dict[id]["canonical_id"])
                            )
            mod_fasta_list.append(record)
        SeqIO.write(mod_fasta_list, mod_fasta_path, "fasta")
class FetchOmaSeq():
    """
    舊方法，抓的比較慢
    get paralogs by uniprot id from OMA, 
    https://omabrowser.org/oma/home/
    """
    def __init__(self):
        pass

    def get_oma_seq(self, uniprot_id, path):
        """
        get paralogs from OMA by uniprot id

        uniprot_id: str, uniprot id
        path: str, path to save fasta file

        return: None
        """
        
        path = Path(path)
        fasta_path = path / "{}.fasta".format(uniprot_id)

        #get orthologs
        orthologs_list = self.__get_orthologs(uniprot_id)

        #writing to fasta
        self.__get_fasta(orthologs_list, fasta_path)

        __uniprot_id_consistance_check(fasta_path, uniprot_id)

    def __get_protein_info_from_entry(self, ortholog_entry):

        try:
            resp = requests.get("https://omabrowser.org/api/protein/{}/".format(ortholog_entry))
            oma_raw = json.loads(resp.text)

            species = oma_raw["species"]["species"]
            # species name is too log, remove some strain info
            species = re.sub("\(.*\)", "", species)

            oma_id = oma_raw["omaid"]
            canonical_id = oma_raw["canonicalid"]
            taxon_id = oma_raw["species"]["taxon_id"]
            sequence = oma_raw["sequence"]

            return {
                "species": species,
                "oma_id": oma_id,
                "canonical_id": canonical_id,
                "taxon_id": taxon_id,
                "sequence": sequence,
            }
        except:
            raise Exception("get single ortholog entry {} from OMA failed".format(ortholog_entry))

    def __get_orthologs(self, uniprot_id):

        try:
            resp = requests.get("https://omabrowser.org/api/group/{}/".format(uniprot_id))
            oma_raw = json.loads(resp.text)

            orthologs_list = []
            t = trange(len(oma_raw["members"]), desc=uniprot_id, leave=True, position=2)

            for i in t:
                orthologs_list.append(self.__get_protein_info_from_entry(oma_raw["members"][i]["entry_nr"]))  
            return orthologs_list
        except:
            raise Exception("get ortholog {} from OMA failed".format(uniprot_id))

    def __get_fasta(self, orthologs_list, path):
        fasta_list = []
        for i in orthologs_list:
            record = SeqRecord(
                Seq(i["sequence"]),
                id=i["oma_id"],
                description="| {} | {} | {}".format(
                    i["species"], i["taxon_id"], i["canonical_id"]
                ),
            )
            fasta_list.append(record)
        SeqIO.write(fasta_list, path, "fasta")


class TaxSeqFilter():
    """
    filter paralogs by taxonomy id, and save as fasta file
    """
    def __init__(self, taxonomy):
        """
        taxonomy: int, taxonomy id from NCBI for filter
                       NCBI: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=9606
        """
        resp = requests.get("https://omabrowser.org/api/taxonomy/{}".format(taxonomy))
        self.taxonomy = taxonomy
        self.taxonomy_list = resp.text

    def taxfilter(self, oma_fasta_path, grouped_fasta_path):
        """
        oma_fasta_path: str, fasta file path for all OMA paralogs
        grouped_fasta_path: str, fasta file path for grouped paralogs
        
        return: None
        """
        # read
        oma_fasta_list = fasta_to_seqlist(oma_fasta_path)

        # filter
        filtered_list = []
        for i in oma_fasta_list:
            tax_id = i.description.split("|")[2].replace(" ", "")
            if tax_id in self.taxonomy_list:
                filtered_list.append(i)

        with open(grouped_fasta_path, "w") as output_handle:
            SeqIO.write(filtered_list, output_handle, "fasta")

    def taxfilter_pipeline(self, oma_path, grouped_tax_path):
        tax_id = Path(str(self.taxonomy))
        species = get_taxid_dict()[self.taxonomy]
        grouped_tax_path.mkdir(exist_ok=True, parents=True)

        fasta_pathlist = list(Path(oma_path).rglob("*.fasta"))
        t = trange(len(fasta_pathlist), desc="", leave=True, position=2)

        for i in t:
            uniprot_id_path = fasta_pathlist[i]
            t.set_description("{} ({}): {}".format(species, tax_id, uniprot_id_path.parts[-1].split(".")[0]))
            self.taxfilter(uniprot_id_path, grouped_tax_path / Path(uniprot_id_path.parts[-1]))
