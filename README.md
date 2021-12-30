# The return of the rings
The return of the rings: evolutionary role of aromatic residues in liquid-liquid phase separation

## Dependencies
* [python 3.9.6](https://www.python.org/)
* [chrome](https://www.google.com/chrome/)
* [Clustal Omega](http://www.clustal.org/omega/)

## Usage
### Clone and build environment
* Clone this repository
    ```
    git clone https://github.com/allmwh/the_return_of_the_rings.git
    ```
* Build environment
    ```
    pip install -r requirements.txt
    ```

After that, just type `jupyter-notebook` in folder, and you can run all of the code
```
jupyter-notebook
``` 
### PONDR order/disorder identification
> [pondr_disorder_ident.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/pondr_disorder_ident.ipynb)   

Get order/disorder identification infos by PONDR 

### Get taxonomy sequence from OMA
> [get_taxonomy_sequence.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/get_taxonomy_sequence.ipynb)   
  
The script does following steps:
* Get paralogs from OMA by uniprot id, downloaded sequences are in `./output/fasta/a_oma`
* Group paralogs by [taxonomy id](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=9606), outputs: `./output/fasta/b_grouped/TAXONOMY_ID`
* Align sequences by Clustal Omega, outputs: `./output/fasta/c_alied/TAXONOMY_ID`
* Filter sequences, outputs: `./output/fasta/d_extre_filtered/TAXONOMY_ID`

### Generate data and figures

> [fig1a.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/fig1a.ipynb)   
> [fig1b.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/fig1b.ipynb)   
> [fig2.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/fig2.ipynb)   

Following scripts instructions to generate data and figures
## Citation
```
@article
TO-DO
```
