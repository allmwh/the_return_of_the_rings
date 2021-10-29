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


class FetchOmaSeq:
    """
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
        try:
            path = Path(path)
            fasta_path = path / "{}.fasta".format(uniprot_id)
            orthologs_list = self.__get_orthologs(uniprot_id)
            self.__get_fasta(orthologs_list, fasta_path)
        except:
            print("{} OMA DOES NOT HAVE THIS UNIPROT_ID RECORD".format(uniprot_id))
            raise Exception("error")

        # some uniprot id in OMA paralogs is not consist with uniprot 
        uniprot_id_oma_fassta = find_human_sequence(fasta_path)["uniprot_id"]
        if uniprot_id != uniprot_id_oma_fassta:
            fasta_path.unlink()
            print("{} IN UNIPROT IS NOT CONSIST WITH OMA's record".format(uniprot_id))
            raise Exception("error")

    def pipeline_get_oma_seq(self, uniprot_ids, path):
        path = Path(path)
        t = trange(len(uniprot_ids), leave=True, position=0)
        for i in t:
            self.get_oma_seq(uniprot_ids[i], path)

    def __get_protein_info_from_entry(self, uniprot_id):
        resp = requests.get("https://omabrowser.org/api/protein/{}/".format(uniprot_id))
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

    def __get_orthologs(self, uniprot_id):
        resp = requests.get("https://omabrowser.org/api/group/{}/".format(uniprot_id))
        oma_raw = json.loads(resp.text)

        orthologs_list = []
        t = trange(len(oma_raw["members"]), desc=uniprot_id, leave=True, position=2)

        for i in t:
            orthologs_list.append(self.__get_protein_info_from_entry(oma_raw["members"][i]["entry_nr"]))  
        return orthologs_list

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


class TaxSeqFilter:
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
