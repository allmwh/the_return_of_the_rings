import re
import numpy as np
from Bio import SeqIO

from function.utilities import fasta_to_seqlist
from function.utilities import find_human_sequence
from function.utilities import get_fasta_seq_info


class FastaExtreFilter:
    """
    filter some homologous, 
    i.e. some sequences in homologous have very long gap in alied fasta file, 
    while other sequences do not have gap, this is "special", and remove it
    """
    def __inti__(self):
        pass

    def fasta_extre_filter(self, input_fasta_path, output_fasta_path, gap_filter_num=20, max_filter_seq=3):
        """
        pepeline： get seq array >> loop(filter by_seq_num >> filter by_gap_num) >> delete extre seq by index
        
        input fasta file, filter special sequence, and output to new fasta

        input_fasta_path: str, fasta file path
        output_fasta_path: str, filtered fasta file path
        gap_filter_num: int, sequence that gap larger than gap_filter_num is possible special sequence
        max_filter_seq: int, sequences number that have same gap simultaneously are special sequence
        """
        #get seq array
        seq_array = self.__seq_2_array(input_fasta_path)
        #get extre seqeunce index
        extre_index = self.__get_extre_index_max_cond(seq_array, gap_filter_num, max_filter_seq)
        #remove sequence by index
        self.__fasta_remove_extre_by_index(input_fasta_path, output_fasta_path, extre_index)

        #error handle for fasta file is empty
        uniprot_id = get_fasta_seq_info(input_fasta_path)['human_uniprot_id']
        if fasta_to_seqlist(output_fasta_path) == list():
            raise Exception("{} filtered fasta is empty".format(uniprot_id))

        return {'before':seq_array.shape[0],
                'after':seq_array.shape[0] - len(extre_index),
                'removed_sequence_index':extre_index}

    def __seq_2_array(self, fasta_path):
        '''
        sequence to numpy array with 1(gap) and 0(residue)

        fasta_path: str, fasta file path

        return: sequence array with 1 and 0
        '''
        seq_info_list = fasta_to_seqlist(fasta_path)
        seq_list = []

        for seq_info in seq_info_list:
            seq = str(seq_info.seq)
            seq = re.sub("[A-Z]", "0", seq)
            seq = re.sub("\-", "1", seq)
            seq = list(seq)
            seq_list.append(seq)
        seq_array = np.array(seq_list, dtype=np.uint8)

        return seq_array

    #### translate annotation TO-DO :) #####
    def __get_extre_index_max_cond(self, seq_array, filter_gap_num, max_filter_seq_num):
        '''
        loop max_filter_seq_num，依序進行get_pot_seq_by_seq_num和get_pot_seq_by_gap_num

        seq_array: np.array, __seq_2_array做的array
        filter_gap_num: int, 多於幾個gap以上是特例
        max_filter_seq_num: int, 少於幾個序列以下是特例

        return: extre_index, 要刪掉的序列index，已經倒序過
        '''
        extre_index = set()

        for seq_num in range(1, max_filter_seq_num + 1):
            extre_array = self.__get_pot_seq_by_seq_num(seq_array, seq_num)
            extre_list = self.__get_pot_seq_by_gap_num(extre_array, filter_gap_num)
            extre_index = extre_index | extre_list
        extre_index = sorted(extre_index, reverse=True)

        return extre_index

    def __get_pot_seq_by_seq_num(self, seq_array, filter_seq_num):
        '''
        column-wise: 看一次出現多少個"特別"在所有序列中，
        一次出現filter_seq_num的序列保留是可能會被刪的，
        array中填回1，填完的array再之後判斷長度

        seq_array: np.array, __seq_2_array做的array
        filter_seq_num: int, 出現"特別"幾次算是特例，要被挑出來
        '''
        extre_array = np.zeros_like(seq_array)
        for seq_index, element in enumerate(seq_array.T):
            #https://li.allmwh.org/8FHWy9
            #參照上面方法，找出unique的，也就是0or1，並得到0or1的index
            vals, inverse, count = np.unique(element, return_inverse=True, return_counts=True)
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

    def __get_pot_seq_by_gap_num(self, extre_array, filter_gap_num):
        '''
        接get_pot_seq_filter_by_seq_num的array，判斷長度，大於filter_gap_num就拿出他的index

        extre_array: np.array, __get_pot_seq_by_seq_num產生的array
        filter_gap_num: 大於多少長度要拿出index
        '''
        problem_list = []
        for row, seq in enumerate(extre_array):
            seq = seq.tolist()
            seq = "".join(str(e) for e in seq)
            gap_num_check = bool(
                re.findall("1{" + re.escape(str(filter_gap_num)) + ",}", seq)
            )
            if gap_num_check:
                problem_list.append(row)

        return set(problem_list)

    def __fasta_remove_extre_by_index(self, input_fasta_path, output_fasta_path, extre_index):
        '''
        輸入input_fasta_path，和想要減掉的序列的index，減掉之後輸出成output_fasta_path

        input_fasta_path: str, fasta路徑
        output_fasta_path: str, fasta路徑
        extre_index: list或可iter的東西，裡面是數字 {1,34,55}
        '''
        seq_info_list = fasta_to_seqlist(input_fasta_path)
        for index in sorted(extre_index, reverse=True): 
            del seq_info_list[index]
        SeqIO.write(seq_info_list, output_fasta_path, "fasta")


class SeqFilter:
    """
    filter order/disorder sequence is shorter than desired length

    od_ident：
        1 disorder
        0 order

        x filtered disorder 
        z filtered order
    """
    def __init__(self):
        pass

    def length_filter_by_od_ident(self, od_ident, disorder_filter_length, order_filter_length):
        '''
        filter order/disorder seqeunce is shorter than order_filter_length/disorder_filter_length

        od_ident: str
        disorder_filter_length: int, continuous disorder sequence which is shorter than disorder_filter_length is replaced with "x" 
        order_filter_length: int, continuous order sequence which is shorter disorder_filter_length is replaced with "z"

        return od_ident
        '''
        #filter_length
        disorder_check = re.finditer("1+", od_ident)
        for i in disorder_check:
            if i:
                start = i.start()
                end = i.end()
                if end - start <= disorder_filter_length:
                    od_ident = od_ident[:start] + "x" * (end - start) + od_ident[end:]
        
        order_check = re.finditer("0+", od_ident)
        for i in order_check:
            if i:
                start = i.start()
                end = i.end()
                if end - start <= order_filter_length:
                    od_ident = od_ident[:start] + "z" * (end - start) + od_ident[end:]

        return od_ident

    def get_seq_from_od_ident(self, od_ident,sequence,od):
        '''
        make new sequence from length filtered od_ident
        
        od_ident: str
        sequence: str, protein seqeunce 
        od: str, ['order', 'disorder'] make order or disorder sequence 
        '''
        new_seq = ''
        for index, element in enumerate(od_ident):
            if od == 'order': # make order seq
                if element == '0':
                    new_seq = new_seq + sequence[index]
                elif element == 'z':
                    new_seq = new_seq + 'z'
                elif element == '1' or element == 'x':
                    new_seq = new_seq + '*'
            
            elif od == 'disorder':# make disorder seq
                if element == '1':
                    new_seq = new_seq + sequence[index]
                elif element == 'x':
                    new_seq = new_seq + 'x'
                elif element == '0' or element == 'z':
                    new_seq = new_seq + '*'
                    
        return new_seq

    def get_od_index(self, od_ident):
        '''
        get order/disorder index

        od_ident: str

        return order/disorder index   
        '''
        order_region = []
        disorder_region = []

        disorder_check = re.finditer("1+", od_ident)
        for i in disorder_check:
            if i:
                start = i.start()
                end = i.end()

                disorder_region.append({"start": start, "end": end})

        order_check = re.finditer("0+", od_ident)
        for i in order_check:
            if i:
                start = i.start()
                end = i.end()

                order_region.append({"start": start, "end": end})

        return {"order_region": order_region, 
                "disorder_region": disorder_region}

    def od_add_alignment(self, path, od_ident):
        '''
        add gap to od_ident to consist with alied sequence length

        path: file, alied fasta 
        od_ident: str

        return: od_ident with gap added
                
        od condition to fill gap
            een: ---(0or1)
            nen0: 0---0
            nen1: 1---1
            nee: (1or0)---
            ned: 1---0 or 0---1
            oeo: (o or empty)---(o or empty)
        '''
        
        # get alied human sequence
        alied_sequence = find_human_sequence(path)["sequence"]
        
        # make new od_ident
        for index, element in enumerate(alied_sequence):
            if element == "-":
                od_ident = od_ident[:index] + "-" + od_ident[index:]


        # both sides are same in center
        # nen1: 1---1
        for _ in range(2):
            nen1_check = re.finditer("1\-+1", od_ident)
            for i in nen1_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])

        # nen0: 0---0
        for _ in range(2):
            nen0_check = re.finditer("0\-+0", od_ident)
            for i in nen0_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])
                    
        # nenx: x---x
        for _ in range(2):
            nenx_check = re.finditer("x\-+x", od_ident)
            for i in nenx_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])
                    
        # nenz: z---z
        for _ in range(2):
            nenz_check = re.finditer("z\-+z", od_ident)
            for i in nenz_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + filed_od * (end - start) + od_ident[end:])


        # start and and to fill 0(order) ,1(disorder), x(disorder filtered) or z(order filtered)
        # een: ---(0, 1, x, z)
        een_check = re.search("^\-+[0,1,x,z]", od_ident)
        if een_check:
            filed_od = een_check.group()[-1]
            start = een_check.start()
            end = een_check.end()
            od_ident = filed_od * (end - start) + od_ident[end:]
            
        #nee: (0, 1, x, z)---
        nee_check = re.search("[0,1,x,z]\-+$", od_ident)
        if nee_check:
            filed_od = nee_check.group()[0]
            start = nee_check.start()
            end = nee_check.end()
            od_ident = od_ident[:start] + filed_od * (end - start)


        # complicated condition
        # nen01: (0, 1)---(0, 1)
        for _ in range(2):
            nen01_check = re.finditer("[0,1]\-+[0,1]", od_ident)
            od_place = "0" 
            for i in nen01_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
                    
        # nenxz: (x, z)---(x, z)
        for _ in range(2):
            nenxz_check = re.finditer("[x,z]\-+[x,z]", od_ident)
            od_place = "z" 
            for i in nenxz_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
                    
        # nen1z: (1, z) --- (z, 1)
        for _ in range(2):
            nen1z_check = re.finditer("1\-+z", od_ident)
            od_place = "z" 
            for i in nen1z_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
        for _ in range(2):
            nen1z_check = re.finditer("z\-+1", od_ident)
            od_place = "z" 
            for i in nen1z_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
                
        # nen0x: (0, x) --- (x, 0)
        for _ in range(2):
            nen0x_check = re.finditer("0\-+x", od_ident)
            od_place = "x" 
            for i in nen0x_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
        for _ in range(2):
            nen0x_check = re.finditer("x\-+0", od_ident)
            od_place = "x" 
            for i in nen0x_check:
                if i:
                    filed_od = i.group()[0]
                    start = i.start()
                    end = i.end()
                    od_ident = (od_ident[:start] + od_place * (end - start) + od_ident[end:])
        
        return od_ident
