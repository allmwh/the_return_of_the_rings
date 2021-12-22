import math
import shutil
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path

from function.utilities import seq_aa_check
from function.utilities import fasta_to_seqlist
from function.utilities import find_human_sequence
from function.utilities import get_uniprot_id_from_fasta

from selenium import webdriver
from selenium.webdriver.common.by import By 
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager


class Rate4Site():
    
    def __init__(self, work_path = '/home/wenlin/tmp'):
        '''
        Calculate conserve score by Rate4Site
        Please install Rate4Site: https://launchpad.net/ubuntu/xenial/+package/rate4site
        Ubuntu/Debian based distro: sudo apt install rate4site

        work_path: str, path to save Rate4Site output file
        '''
        
        work_path = Path(work_path)
        self.work_path = work_path
        
        self.res_path = work_path / "tmp.res"   
    
    def get_conserve_score(self, fasta_path):
        
        #run rate4site
        self.__run_rate4site(fasta_path, self.res_path)
        
        #parse_res_file
        conserve_score = self.__parse_res_file(self.res_path)

        #delete r4s output 
        shutil.rmtree(str(self.work_path))
        self.work_path.mkdir(parents=True, exist_ok=True)

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

    def __init__(self):
        '''
        Please consider use standalone version ConserveStandalone() to have better speed instead of web crawler 

        Claculate conserve score by web
        https://compbio.cs.princeton.edu/conservation/score.html
        '''
        options = webdriver.ChromeOptions()
        
        options.add_argument('--headless')
        options.add_argument('--no-sandbox')
        options.add_argument("--disable-extensions")
        
        self.driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)

    def get_conserve_score(self, fasta_path, method, windows_size=0, sequence_weight=False):
        '''
        windows_size: int
        sequence_weighting: bool
        method: ['jsd', 'shannon_entropy', 'property_entropy', 'relative_entropy', 'sum_of_pairs']
        fasta_path = alied過的fasta檔名
        
        傳入單一個檔案，輸出分數list
        '''
        base_url = 'https://compbio.cs.princeton.edu/conservation/score.html'
        self.driver.get(base_url)
        fasta_path = fasta_path.resolve()

        #windows size
        windows_size_buttom = self.driver.find_element(By.XPATH,'/html/body/div/div[2]/form/input[1]')
        windows_size_buttom.clear()
        windows_size_buttom.send_keys(windows_size) #要先清除才能送keys
        
        #Sequence Weighting
        sequence_weighting_bottom = self.driver.find_element(By.XPATH,'/html/body/div[1]/div[2]/form/input[2]')
        sequence_weighting_bottom.click() #clean default click
        if sequence_weight == True:
            sequence_weighting_bottom.click()
        elif sequence_weight == False:
            pass
            
        #choose method
        self.__choose_method(method)
        
        #send data
        self.driver.find_element(By.XPATH, '/html/body/div/div[2]/form/p[5]/input').send_keys(str(fasta_path))
        self.driver.find_element(By.XPATH, '/html/body/div/div[2]/form/p[7]/input').click()
        
        #wait
        webdriver.support.ui.WebDriverWait(self.driver, 2).until(EC.presence_of_element_located((By.XPATH, "/html/body/a")))
        
        #save result to list
        result = self.driver.find_element(By.XPATH,'/html/body/pre')
        score = []
        for i in result.text.split('\n')[2:]:
            i = i.split(' ')
            value = float(i[1])
            score.append(value)
        
        return score

    def __choose_method(self, method):
        method_select = Select(self.driver.find_element(By.XPATH,'/html/body/div[1]/div[2]/form/select[1]'))
        method_dict= {'jsd':lambda: method_select.select_by_visible_text('JS Divergence'),
                      'shannon_entropy':lambda: method_select.select_by_visible_text('Shannon entropy'),
                      'property_entropy':lambda: method_select.select_by_visible_text('Property entropy'),
                      'vn_entropy':lambda: method_select.select_by_visible_text('VN entropy'),
                      'relative_entropy':lambda: method_select.select_by_visible_text('Relative entropy'),
                      'sum_of_pairs':lambda: method_select.select_by_visible_text('Sum of Pairs')}
        method_dict[method]()

    def close_webdriver(self):
        return self.driver.close()      


class ConserveStandalone():
    
    def __init__(self):
        '''
        Claculate conserve score by standalone code by 
        Capra JA and Singh M. Predicting functionally important residues from sequence conservation. 
        Bioinformatics, 23(15):1875-82, 2007. [Bioinformatics]
        https://compbio.cs.princeton.edu/conservation/score.html
        '''
        self.__PSEUDOCOUNT = 0.0000001
        self.__amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
        self.__blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
        aa_to_index = {}
        for i, aa in enumerate(self.__amino_acids):
            aa_to_index[aa] = i
        self.__aa_to_index = aa_to_index

    def get_conserve_score(self, fasta_path, method, windows_size=0, sequence_weight=False):
        '''
        windows_size: int
        sequence_weighting: bool
        method: ['jsd', 'shannon_entropy']
        fasta_path = alied過的fasta檔名
        
        傳入單一個檔案，輸出分數list
        '''
        pass

        seqs_strip = self.__get_seqs_strip(fasta_path)
        gap_penalty=1

        #sequence_weight
        if sequence_weight:
            weight = self.__calculate_sequence_weights(seqs_strip)
        else:
            weight = [1] #no sequence weight

        score_list = []
        for col in seqs_strip:
            
            #method choose
            if method == 'shannon_entropy':
                score = self.__shannon_entropy(col=col, 
                                            seq_weights=weight, 
                                            gap_penalty=gap_penalty)
            elif method == 'jsd':
                score = self.__js_divergence(col=col, 
                                            bg_distr=self.__blosum_background_distr, 
                                            seq_weights=weight,
                                            gap_penalty=gap_penalty)

            score_list.append(score)

        #windows size
        if windows_size > 0:
            score_list = self.__window_score(score_list, windows_size)

        return score_list

    def __get_seqs_strip(self, fasta_path):

        fasta_path = Path(fasta_path)
        fasta_list = fasta_to_seqlist(fasta_path)

        all_seq = [list(seq_aa_check(str(seqrecord.seq))) for seqrecord in fasta_list]
        seqs_strip = np.array(all_seq, dtype=object).T
        
        return seqs_strip

    def __weighted_freq_count_pseudocount(self, col, seq_weights):
        """ Return the weighted frequency count for a column--with pseudocount."""

        # if the weights do not match, use equal weight
        if len(seq_weights) != len(col):
            seq_weights = [1.] * len(col)

        aa_num = 0
        freq_counts = len(self.__amino_acids)*[self.__PSEUDOCOUNT] # in order defined by amino_acids

        for aa in self.__amino_acids:
            for j in range(len(col)):
                if col[j] == aa:
                    freq_counts[aa_num] += 1 * seq_weights[j]

            aa_num += 1

        for j in range(len(freq_counts)):
            freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(self.__amino_acids) * self.__PSEUDOCOUNT)

        return freq_counts

    def __calculate_sequence_weights(self, seqs_strip):
        """ Calculate the sequence weights using the Henikoff '94 method
        for the given msa. """

        seqs_strip = seqs_strip.T #return back to msa form

        seq_weights = [0.] * len(seqs_strip)
        for i in range(len(seqs_strip[0])):
                freq_counts = [0] * len(self.__amino_acids)

                col = []
                for j in range(len(seqs_strip)):
                    if seqs_strip[j][i] != '-': # ignore gaps
                        freq_counts[self.__aa_to_index[seqs_strip[j][i]]] += 1

                num_observed_types = 0
                for j in range(len(freq_counts)):
                    if freq_counts[j] > 0: num_observed_types +=1

                for j in range(len(seqs_strip)):
                    d = freq_counts[self.__aa_to_index[seqs_strip[j][i]]] * num_observed_types
                    if d > 0:
                        seq_weights[j] += 1. / d

        for w in range(len(seq_weights)):
                seq_weights[w] /= len(seqs_strip[0])

        return seq_weights

    def __weighted_gap_penalty(self, col, seq_weights):
        """ Calculate the simple gap penalty multiplier for the column. If the 
        sequences are weighted, the gaps, when penalized, are weighted 
        accordingly. """

        # if the weights do not match, use equal weight
        if len(seq_weights) != len(col):
            seq_weights = [1.] * len(col)
        
        gap_sum = 0.
        for i in range(len(col)):
            if col[i] == '-':
                gap_sum += seq_weights[i]

        return 1 - (gap_sum / sum(seq_weights))
        
    def __shannon_entropy(self, col, seq_weights, gap_penalty=1):

        """Calculates the Shannon entropy of the column col. sim_matrix  and 
        bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The 
        entropy will be between zero and one because of its base. See p.13 of 
        Valdar 02 for details. The information score 1 - h is returned for the sake 
        of consistency with other scores."""

        fc = self.__weighted_freq_count_pseudocount(col, seq_weights)

        h = 0. 
        for i in range(len(fc)):
            if fc[i] != 0:
                h += fc[i] * math.log(fc[i])

        h /= math.log(min(len(fc), len(col)))

        inf_score = 1 - (-1 * h)

        if gap_penalty == 1: 
            return round((inf_score * self.__weighted_gap_penalty(col, seq_weights)),5)
        else: 
            return round(inf_score,5)

    def __js_divergence(self, col, bg_distr, seq_weights, gap_penalty=1):
        """ Return the Jensen-Shannon Divergence for the column with the background
        distribution bg_distr. sim_matrix is ignored. JSD is the default method.
        """
        
        distr = bg_distr[:]

        fc = self.__weighted_freq_count_pseudocount(col, seq_weights)

        # if background distrubtion lacks a gap count, remove fc gap count
        if len(distr) == 20: 
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in range(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc

        if len(fc) != len(distr): return -1

        # make r distriubtion
        r = []
        for i in range(len(fc)):
            r.append(.5 * fc[i] + .5 * distr[i])

        d = 0.
        for i in range(len(fc)):
            if r[i] != 0.0:
                if fc[i] == 0.0:
                    d += distr[i] * math.log(distr[i]/r[i], 2)
                elif distr[i] == 0.0:
                    d += fc[i] * math.log(fc[i]/r[i], 2) 
                else:
                    d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)

        # d /= 2 * math.log(len(fc))
        d /= 2

        if gap_penalty == 1: 
            return round((d * self.__weighted_gap_penalty(col, seq_weights)),5) 
        else:  
            return round(d,5) 
        
    def __window_score(self, scores, window_len, lam=.5):
        """ This function takes a list of scores and a length and transforms them 
        so that each position is a weighted average of the surrounding positions. 
        Positions with scores less than zero are not changed and are ignored in the 
        calculation. Here window_len is interpreted to mean window_len residues on 
        either side of the current residue. """

        w_scores = scores[:]

        for i in range(window_len, len(scores) - window_len):
            if scores[i] < 0: 
                continue

            sum = 0.
            num_terms = 0.
            for j in range(i - window_len, i + window_len + 1):
                if i != j and scores[j] >= 0:
                    num_terms += 1
                    sum += scores[j]

            if num_terms > 0:
                w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]

        return [round(i,5) for i in w_scores]


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
        human_sequence = find_human_sequence(fasta_path)["remove_gap_sequence"]

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

        
        #length check
        uniprot_id = get_uniprot_id_from_fasta(fasta_path)
        if not (len(human_sequence) == len(conserve_score) == len(od_ident)):
            raise Exception('''{} human_sequence, conserve_score, od_ident length check, these three length must be same,
            possible not same reason: 
            1.OMA database with same uniprot_id, while sequence is different
            2.human_sequence does not remove "-", because od_ident and conserve_score is calculated by remove gap
            '''.format(uniprot_id))

        for aa, score, od in zip(human_sequence, conserve_score, od_ident):
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
                disorder_score_dict[key] = round((disorder_score_dict[key] / disorder_content_dict[key]), 5)
            if order_content_dict[key] == 0:
                order_score_dict[key] = 0
            else:
                order_score_dict[key] = round((order_score_dict[key] / order_content_dict[key]), 5)

        # content sum
        order_content_dict["total"] = sum(order_content_dict.values())
        disorder_content_dict["total"] = sum(disorder_content_dict.values())

        return {
            "conserve": {"order": order_score_dict, "disorder": disorder_score_dict},
            "content": {"order": order_content_dict, "disorder": disorder_content_dict},
        }
