from pathlib import Path
from tqdm.notebook import trange
from Bio.Align.Applications import ClustalOmegaCommandline


class Alignment:
    def __init__(self):
        pass

    def alignment_single(self, input_fasta, output_fasta):
        """
        single fasta file alignment

        input_fasta: fasta path for input, Ex: /home/wenlin/d/rbp/oma/base/Q13148.fasta
        output_fasta: fasta path for output, Ex: /home/wenlin/d/rbp/oma/base/new.fasta
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

        input_path: path contains many fasta files, Ex: /home/wenlin/d/rbp/oma/base/
        output_path: path contains many fasta files for output, Ex: /home/wenlin/d/rbp/oma/new/
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
