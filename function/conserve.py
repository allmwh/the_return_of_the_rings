import math
import numpy as np
from pathlib import Path

from function.utilities import fasta_to_seqlist
from function.utilities import find_human_sequence
from function.utilities import get_only_human_score
from function.utilities import get_uniprot_id_from_fasta

from selenium import webdriver
from selenium.webdriver.common.by import By 
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

from tqdm.notebook import trange



import subprocess
import pandas as pd
from pathlib import Path
from function.utilities import find_human_sequence

class Rate4Site():
    def __init__(self, work_path = '/home/wenlin/tmp'):
        
        work_path = Path(work_path)
        self.work_path = work_path
        
        self.res_path = work_path / "tmp.res"
        
    
    def get_conserve_score(self, fasta_path):
        
        #run rate4site
        self.__run_rate4site(fasta_path, self.res_path)
        
        #parse_res_file
        conserve_score = self.__parse_res_file(self.res_path)
        
        return conserve_score
        
    def __run_rate4site(self, fasta_path,res_path):
        

        fasta_path = str(fasta_path)
        res_path = str(res_path)

        #get ref seq name
        seqence_name = find_human_sequence(fasta_path)['sequence_name']
        
        
        #-s input fasta
        #-o output score 
        #-a reference sequence name, human sequence's name in fasta
        try:
            process = subprocess.Popen(['rate4site', '-s', fasta_path, '-o', res_path, '-a', seqence_name],
                                        cwd = self.work_path,
                                        stdout = subprocess.PIPE, 
                                        stderr = subprocess.PIPE)
            stdout, stderr = process.communicate()
        except:
            raise FileNotFoundError("Rate4Site may not be installed, please install Rate4Site")

        
    def __parse_res_file(self, res_path):
        df = pd.read_csv(res_path,skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,12], 
                 delim_whitespace=True, 
                 skipfooter=2, 
                 header=None,
                 engine='python')
        return df[2].to_list()


class ConserveByWeb():
    '''
    從這網站丟序列抓保守分數，有很多不同的算法
    https://compbio.cs.princeton.edu/conservation/score.html
    '''
    def __init__(self):
        options = webdriver.ChromeOptions()
        
        options.add_argument('--headless')
        options.add_argument('--no-sandbox')
        options.add_argument("--disable-extensions")
        # options.add_argument("--disable-gpu")
        # options.add_argument('--disable-dev-shm-usage')
        
        self.driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)

    def close(self):
        return self.driver.close()

    def choose_method(self,method):
        method_select = Select(self.driver.find_element(By.XPATH,'/html/body/div[1]/div[2]/form/select[1]'))
        method_dict= {'jsd':lambda: method_select.select_by_visible_text('JS Divergence'),
                      'shannon_entropy':lambda: method_select.select_by_visible_text('Shannon entropy'),
                      'property_entropy':lambda: method_select.select_by_visible_text('Property entropy'),
                      'vn_entropy':lambda: method_select.select_by_visible_text('VN entropy'),
                      'relative_entropy':lambda: method_select.select_by_visible_text('Relative entropy'),
                      'sum_of_pairs':lambda: method_select.select_by_visible_text('Sum of Pairs')}
        method_dict[method]()
        
        
    def calculate(self,windows_size,sequence_weighting,method,alied_fasta):
        '''
        windows_size: int
        sequence_weighting: bool
        method: ['jsd', 'shannon_entropy', 'property_entropy', 'relative_entropy', 'sum_of_pairs']
        alied_fasta = alied過的fasta檔名
        
        傳入單一個檔案，輸出分數list
        '''
        base_url = 'https://compbio.cs.princeton.edu/conservation/score.html'
        self.driver.get(base_url)

        #windows size
        windows_size_buttom = self.driver.find_element(By.XPATH,'/html/body/div/div[2]/form/input[1]')
        windows_size_buttom.clear()
        windows_size_buttom.send_keys(windows_size) #要先清除才能送keys
        
        #Sequence Weighting
        sequence_weighting_bottom = self.driver.find_element(By.XPATH,'/html/body/div[1]/div[2]/form/input[2]')
        sequence_weighting_bottom.click() #clean default click
        if sequence_weighting == True:
            sequence_weighting_bottom.click()
        elif sequence_weighting == False:
            pass
            
        #choose method
        self.choose_method(method)
        
        #send data
        self.driver.find_element(By.XPATH, '/html/body/div/div[2]/form/p[5]/input').send_keys(str(alied_fasta))
        self.driver.find_element(By.XPATH, '/html/body/div/div[2]/form/p[7]/input').click()
        
        #wait
        webdriver.support.ui.WebDriverWait(self.driver, 2).until(EC.presence_of_element_located((By.XPATH, "/html/body/a")))
        
        #save result to list
        result = self.driver.find_element(By.XPATH,'/html/body/pre')

        score = []
        for i in result.text.split('\n')[2:]:

            i = i.split(' ')

            index = i[0]
            value = float(i[1])
            seq = i[2]

            score.append(value)
        
        return score

class Entropy():
    def __init__(self):
        """
        Capra JA and Singh M. Predicting functionally important residues from sequence conservation.
        Bioinformatics, 23(15):1875-82, 2007. [Bioinformatics]
        
        https://compbio.cs.princeton.edu/conservation/
        """
        pass

    def get_conserve_score(self, fasta_path):
        fasta_path = Path(fasta_path)
        fasta_list = fasta_to_seqlist(fasta_path)

        all_seq = [list(str(seqrecord.seq)) for seqrecord in fasta_list]
        all_strip = np.array(all_seq, dtype=object).T
        entropy_list = []

        for strip in all_strip:
            strip = "".join(strip.tolist())
            entropy = self.perstrip_entropy(strip)
            entropy_list.append(entropy)
        return entropy_list

    def perstrip_entropy(self, col, seq_weights=[1], gap_penalty=1):
        """Calculates the Shannon entropy of the column col. sim_matrix  and
        bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The
        entropy will be between zero and one because of its base. See p.13 of
        Valdar 02 for details. The information score 1 - h is returned for the sake
        of consistency with other scores."""
    
        fc = self.__weighted_freq_count_pseudocount(col, seq_weights)
        h = 0.0
        for i in range(len(fc)):
            if fc[i] != 0:
                h += fc[i] * math.log(fc[i])

        # h /= math.log(len(fc))
        h /= math.log(min(len(fc), len(col)))
        inf_score = 1 - (-1 * h)
        if set(col) == {"-"}:
            return np.nan

        if gap_penalty == 1:
            return round((inf_score * self.__weighted_gap_penalty(col, seq_weights)), 3)
        else:
            return round(inf_score, 3)

    def __weighted_freq_count_pseudocount(self, col, seq_weights=[1], pc_amount=0.0000001):
        """Return the weighted frequency count for a column--with pseudocount."""
        amino_acids = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-",]

        # if the weights do not match, use equal weight
        if len(seq_weights) != len(col):
            seq_weights = [1.0] * len(col)
        aa_num = 0
        freq_counts = len(amino_acids) * [pc_amount]  # in order defined by amino_acids
        for aa in amino_acids:
            for j in range(len(col)):
                if col[j] == aa:
                    freq_counts[aa_num] += 1 * seq_weights[j]
            aa_num += 1
        for j in range(len(freq_counts)):
            freq_counts[j] = freq_counts[j] / (
                sum(seq_weights) + len(amino_acids) * pc_amount
            )
        return freq_counts

    def __weighted_gap_penalty(self, col, seq_weights=[1]):
        """Calculate the simple gap penalty multiplier for the column. If the
        sequences are weighted, the gaps, when penalized, are weighted
        accordingly."""
        # if the weights do not match, use equal weight
        if len(seq_weights) != len(col):
            seq_weights = [1.0] * len(col)
        gap_sum = 0.0
        for i in range(len(col)):
            if col[i] == "-":
                gap_sum += seq_weights[i]
        return 1 - (gap_sum / sum(seq_weights))


class ConservePeraa():
    """
    calculate conservation lever per amino acid
    """

    def __init__(self):
        pass



    def get_aa_info(self, fasta_path, conserve_score, od_ident):
        '''
        fasta_path: alied 好的fasta
        conserve_score: conserve score
        od_ident: 切掉非human的od_ident
        '''
        # human sequence
        human_alied_sequence = find_human_sequence(fasta_path)["sequence"]

        order_content_dict = {
            "A": 0, "C": 0, "D": 0,
            "E": 0, "F": 0, "G": 0,
            "H": 0, "I": 0, "K": 0,
            "L": 0, "M": 0, "N": 0,
            "P": 0, "Q": 0, "R": 0,
            "S": 0, "T": 0, "V": 0,
            "W": 0, "Y": 0, "-": 0,
            "total": 0,
            }

        disorder_content_dict = {
            "A": 0, "C": 0, "D": 0,
            "E": 0, "F": 0, "G": 0,
            "H": 0, "I": 0, "K": 0,
            "L": 0, "M": 0, "N": 0,
            "P": 0, "Q": 0, "R": 0,
            "S": 0, "T": 0, "V": 0,
            "W": 0, "Y": 0, "-": 0,
            "total": 0,
            }

        order_score_dict = {
            "A": 0, "C": 0, "D": 0,
            "E": 0, "F": 0, "G": 0,
            "H": 0, "I": 0, "K": 0,
            "L": 0, "M": 0, "N": 0,
            "P": 0, "Q": 0, "R": 0,
            "S": 0, "T": 0, "V": 0,
            "W": 0, "Y": 0, "-": 0,
            }

        disorder_score_dict = {
            "A": 0, "C": 0, "D": 0,
            "E": 0, "F": 0, "G": 0,
            "H": 0, "I": 0, "K": 0,
            "L": 0, "M": 0, "N": 0,
            "P": 0, "Q": 0, "R": 0,
            "S": 0, "T": 0, "V": 0,
            "W": 0, "Y": 0, "-": 0,
            }

        # ERROR handle due to the lehgth of sequence from OMA and uniprot are different, while uniprot_id is the same
        uniprot_id = get_uniprot_id_from_fasta(fasta_path)
        if len(conserve_score) != len(od_ident):
            print("{} ENTROPY LENGTH IS NOT EQUAL WITH OD_IODENT".format(uniprot_id))
            raise Exception("error")

        for aa, score, od in zip(human_alied_sequence, conserve_score, od_ident):
            if np.isnan(score):  
                continue
            if od == "0":  # order
                order_score_dict[aa] = order_score_dict[aa] + score
                order_content_dict[aa] = order_content_dict[aa] + 1
            elif od == "1":  # disorder
                disorder_score_dict[aa] = disorder_score_dict[aa] + score
                disorder_content_dict[aa] = disorder_content_dict[aa] + 1

        # del "-"
        del order_content_dict["-"]
        del disorder_content_dict["-"]
        del order_score_dict["-"]
        del disorder_score_dict["-"]

        # mean
        for key, value in disorder_score_dict.items():
            if disorder_content_dict[key] == 0:
                disorder_score_dict[key] = 0
            else:
                disorder_score_dict[key] = round((disorder_score_dict[key] / disorder_content_dict[key]), 3)
            if order_content_dict[key] == 0:
                order_score_dict[key] = 0
            else:
                order_score_dict[key] = round((order_score_dict[key] / order_content_dict[key]), 3)

        # content sum
        order_content_dict["total"] = sum(order_content_dict.values())
        disorder_content_dict["total"] = sum(disorder_content_dict.values())

        return {
            "conserve": {"order": order_score_dict, "disorder": disorder_score_dict},
            "content": {"order": order_content_dict, "disorder": disorder_content_dict},
        }
