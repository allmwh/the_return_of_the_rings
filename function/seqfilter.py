import re
import numpy as np
from Bio import SeqIO
from itertools import groupby

from function.utilities import find_human_sequence
from function.utilities import fasta_to_seqlist
from function.utilities import get_uniprot_id_from_fasta




class FastaExtreFilter():
    def __inti__(self):
        pass
    
    def fasta_extre_filter(self,input_fasta_path,output_fasta_path,gap_filter_num=20,max_filter_seq=3):
        '''
        pepeline接起來：拿到seq array >> 迴圈(by_seq_num過濾 >> by_gap_num過濾) >> 刪掉有問題的
        '''
        #拿到array
        seq_array = self.__seq_2_array(input_fasta_path)

        #迴圈拿index
        extre_index = self.__get_extre_index_max_cond(seq_array, gap_filter_num, max_filter_seq)

        #刪掉有問題的
        self.__fasta_remove_extre_by_index(input_fasta_path,output_fasta_path,extre_index)

        #ERROR 篩到的序列都沒了
        uniprot_id = get_uniprot_id_from_fasta(input_fasta_path)
        if fasta_to_seqlist(output_fasta_path) == list():
            print("{} FILTERED FASTA IS EMPTY".format(uniprot_id))
            raise Exception("error")

        return extre_index

    def __seq_2_array(self,fasta_path):
        '''
        把序列轉換成array，1是gap，0是序列

        fasta_path: str, fasta路徑

        return: 序列裝在array裡
        '''
        seq_info_list = fasta_to_seqlist(fasta_path)
        seq_list = []
        for seq_info in seq_info_list:
            seq = str(seq_info.seq)
            seq = re.sub("[A-Z]","0",seq)
            seq = re.sub("\-","1",seq)
            seq = list(seq)
            seq_list.append(seq)
        seq_array = np.array(seq_list,dtype=np.uint8)
        return seq_array

    def __get_extre_index_max_cond(self,seq_array, filter_gap_num, max_filter_seq_num):
        '''
        loop max_filter_seq_num，依序進行get_pot_seq_by_seq_num和get_pot_seq_by_gap_num

        seq_array: np.array, __seq_2_array做的array
        filter_gap_num: int, 多於幾個gap以上是特例
        max_filter_seq_num: int, 少於幾個序列以下是特例

        return: extre_index, 要刪掉的序列index，已經倒序過
        '''
        extre_index = set()
        for seq_num in range(1,max_filter_seq_num+1):
            extre_array = self.__get_pot_seq_by_seq_num(seq_array,seq_num)
            extre_list = self.__get_pot_seq_by_gap_num(extre_array,filter_gap_num)
            extre_index = extre_index | extre_list
        extre_index = sorted(extre_index,reverse=True)
        return extre_index
        
    def __get_pot_seq_by_seq_num(self,seq_array,filter_seq_num):
        '''
        column-wise: 看一次出現多少個"特別"在所有序列中，
        一次出現filter_seq_num的序列保留是可能會被刪的，
        array中填回1，填完的array再之後判斷長度

        seq_array: np.array, __seq_2_array做的array
        filter_seq_num: int, 出現"特別"幾次算是特例，要被挑出來
        '''
        extre_array = np.zeros_like(seq_array)
        for seq_index,element in enumerate(seq_array.T):
            #https://li.allmwh.org/8FHWy9
            #參照上面方法，找出unique的，也就是0or1，並得到0or1的index
            vals,inverse,count = np.unique(element, return_inverse=True,return_counts=True)
            idx_vals_repeated = np.where(count >= 1)[0]
            vals_repeated = vals[idx_vals_repeated]
            rows, cols = np.where(inverse == idx_vals_repeated[:, np.newaxis])
            _, inverse_rows = np.unique(rows, return_index=True)
            res = np.split(cols, inverse_rows[1:])

            #res是0 or 1分別的index，先看數量有無符合2，代表不是全0or全1
            #再看符合filter_seq_num數量的是要的index
            #可疑的index填入extre_array中
            if len(res) == 2:
                for seq_nums in res:
                    if len(seq_nums) == filter_seq_num:
                        for seq_num in seq_nums:
                            extre_array[seq_num][seq_index] = 1
        return extre_array
    
    def __get_pot_seq_by_gap_num(self,extre_array,filter_gap_num):
        '''
        接get_pot_seq_filter_by_seq_num的array，判斷長度，大於filter_gap_num就拿出他的index

        extre_array: np.array, __get_pot_seq_by_seq_num產生的array
        filter_gap_num: 大於多少長度要拿出index
        '''
        problem_list = []
        for row,seq in enumerate(extre_array):
            seq = seq.tolist()
            seq = ''.join(str(e) for e in seq)
            gap_num_check = bool(re.findall("1{"+re.escape(str(filter_gap_num))+",}",seq))
            if gap_num_check:
                problem_list.append(row)
        return set(problem_list)


    def __fasta_remove_extre_by_index(self,input_fasta_path,output_fasta_path,extre_index):
        '''
        輸入input_fasta_path，和想要減掉的序列的index，減掉之後輸出成output_fasta_path

        input_fasta_path: str, fasta路徑
        output_fasta_path: str, fasta路徑
        extre_index: list或可iter的東西，裡面是數字 {1,34,55}
        '''
        seq_info_list = fasta_to_seqlist(input_fasta_path)
        for index in sorted(extre_index,reverse=True): #排序從後面刪過來
            del seq_info_list[index] 
        SeqIO.write(seq_info_list, output_fasta_path, "fasta")

class SeqFilter():
    def __init__(self):
        pass
    
    def length_filter_by_od_ident(self,od_ident,disorder_filter_length,order_filter_length):
        #filter_length
        disorder_check = re.finditer("1+",od_ident)
        for i in disorder_check:
            if i:
                start = i.start()
                end = i.end()
                if end-start < disorder_filter_length:
                    od_ident = od_ident[:start] + 'x'*(end-start) + od_ident[end:]
                    
        order_check = re.finditer("0+",od_ident)
        for i in order_check:
            if i:
                start = i.start()
                end = i.end()
                if end-start < order_filter_length:
                    od_ident = od_ident[:start] + 'x'*(end-start) + od_ident[end:]

        return od_ident

    def get_od_index(self,od_ident):
        #get order_disorder index
        order_region = []
        disorder_region = []

        disorder_check = re.finditer("1+",od_ident)
        for i in disorder_check:
            if i:
                start = i.start()
                end = i.end()
                
                disorder_region.append({"start":start,
                                        "end":end})
                
        order_check = re.finditer("0+",od_ident)
        for i in order_check:
            if i:
                start = i.start()
                end = i.end()
                
                order_region.append({"start":start,
                                        "end":end})

        return {"order_region":order_region,
                "disorder_region":disorder_region}
    
    def od_add_alignment(self, path, od_ident):
        '''
        因為od_ident跟alignemt過的長度會不一樣，所以要用alied_seqeunce做新的od_ident

        path: file, alied好的fasta
        od_ident: str, 標示為0or1的disorder判斷句

        return: od_ident,新的o/d判斷，跟alied過的一樣長囉
                
        後面od判斷式
            een: ---(0or1)
            nen0: 0---0
            nen1: 1---1
            nee: (1or0)---
            ned: 1---0 or 0---1
            oeo: (o or empty)---(o or empty)

            再中間的可能會遇到1---1---1這樣的，要補兩次才能把左右邊都補起來
        '''
        
        # 找人的功能割到utilities.find_human_sequence()
        alied_sequence = find_human_sequence(path)["sequence"]
        
        # 做新的od_ident
        for index,element in enumerate(alied_sequence):
            if element == '-':
                od_ident = od_ident[:index] + '-' + od_ident[index:]
        
        
        #一堆判斷
        #een: ---(0or1)
        een_check = re.search('^\-+[0,1,o]',od_ident)
        if een_check:
            filed_od = een_check.group()[-1]
            start = een_check.start()
            end = een_check.end()
            od_ident = filed_od*(end-start) + od_ident[end:]


        #nen0: 0---0
        for _ in range(2):
            nen0_check = re.finditer('0\-+0',od_ident)
            for i in nen0_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = od_ident[:start] + filed_od*(end-start) + od_ident[end:]


        #nen1: 1---1
        for _ in range(2):       
            nen1_check = re.finditer('1\-+1',od_ident)
            for i in nen1_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = od_ident[:start] + filed_od*(end-start) + od_ident[end:]


        #nee: (1or0)---
        nee_check = re.search('[0,1,o]\-+$',od_ident)
        if nee_check:
            filed_od = nee_check.group()[0]
            start = nee_check.start()
            end = nee_check.end()
            od_ident = od_ident[:start] + filed_od*(end-start)


        #ned: 1---0 or 0---1
        for _ in range(2): 
            ned_check = re.finditer('[0,1]\-+[0,1]',od_ident)
            od_place = '0' #兩邊不一樣，先放成order序列(0)
            for i in ned_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = od_ident[:start] + od_place*(end-start) + od_ident[end:]


        #oeo: (o or empty)---(o or empty), fill "o"
        for _ in range(2):
            oeo_check = re.finditer('o?\-+o?',od_ident)
            for i in oeo_check:
                if i:
                    start = i.start()
                    end = i.end()
                    od_ident = od_ident[:start] + "o"*(end-start) + od_ident[end:]

        return od_ident
       