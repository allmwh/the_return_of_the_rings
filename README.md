# The return of the rings
The return of the rings: evolutionary role of aromatic residues in liquid-liquid phase separation

## Dependencies
* [python 3+](https://www.python.org/)
* [chrome](https://www.google.com/chrome/) and [chromedriver](https://chromedriver.chromium.org/)
* [Clustal Omega v1.2.4](http://www.clustal.org/omega/)

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
> [data_pondr_disorder_ident.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/data_pondr_disorder_ident.ipynb)   

Please change `chromedriver_path` to your chrome driver as shown below, and follow instructions in script file 
```
#####CHANGE HERE#####
algorithm = 'VSL2' #PONDR algorithm for use (‘VLXT’, ‘XL1_XT’, ‘CAN_XT’, ‘VL3-BA’, 'VSL2')
chromedriver_path = '/home/wenlin/d/custom_command/chromedriver' #chrome driver path
#####CHANGE HERE#####
```
### Get taxonomy sequence from OMA
> [data_get_taxonomy_seq.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/data_get_taxonomy_seq.ipynb)   
  
This script does following steps:
* Get paralogs from OMA by uniprot id, downloaded sequences are in `./output/fasta/a_oma`
* Group paralogs by [taxonomy id](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=9606), outputs: `./output/fasta/b_grouped/TAXONOMY_ID`
* Align sequences by Clustal Omega, outputs: `./output/fasta/c_alied/TAXONOMY_ID`
* Filter special sequences, outputs: `./output/fasta/d_extre_filtered/TAXONOMY_ID`

### Generate data and figures

> [fig1_a.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/fig1_a.ipynb)   
> [fig1_b.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/fig1_b.ipynb)   
> [fig2.ipynb](https://github.com/allmwh/the_return_of_the_rings/blob/main/fig2.ipynb)   

Following scripts instructions to generate data and figures
## Citation
```
@article
TO-DO
```
