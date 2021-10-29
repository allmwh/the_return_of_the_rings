from pathlib import Path
from tqdm.notebook import trange
from Bio.Align.Applications import ClustalOmegaCommandline


class Alignment:
    """
    alignment func by Biopython ClustalOmegaCommandline
    need to install ClustalOmega: http://www.clustal.org/omega/
    """
    def __init__(self):
        pass

    def alignment_single(self, input_fasta, output_fasta):
        """
        single fasta file alignment

        input_fasta: str, fasta file path for alignment
        output_fasta: str, fasta file path for alied fasta

        return: None
        """
        input_fasta = Path(input_fasta)
        output_fasta = Path(output_fasta)

        try:
            clustalomega_cline = ClustalOmegaCommandline(
                infile=input_fasta,
                outfile=output_fasta,
                verbose=True,
                auto=True,
                force=True,
            )
            clustalomega_cline()
        except:
            print("ALIGNMENT {} failed".format(str(input_fasta)))
            raise Exception("error")

    def alignment_path(self, input_path, output_path):
        """
        alignment all file under input_oath, and output to output_path
        """
        input_path = Path(input_path)
        input_path.mkdir(parents=True, exist_ok=True)
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)

        fasta_pathlist = list(Path(input_path).rglob("*.fasta"))

        t = trange(len(fasta_pathlist), desc="", leave=True, position=0)
        for i in t:
            fasta_path = fasta_pathlist[i].parts[-1]
            t.set_description(fasta_path)

            input_fasta = input_path / fasta_path
            output_fasta = output_path / fasta_path

            self.alignment_single(input_fasta, output_fasta)
