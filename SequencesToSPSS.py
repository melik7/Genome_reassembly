import random
import pprint
from typing import Dict, List, Tuple
from tools_karkkainen_sanders import simple_kark_sort

import time
from timer import Timer 
from progress_bar import update_progress 
import argparse

# We define all the arguments 

parser = argparse.ArgumentParser(description='Build the spss from a read file and calculate FN and FP compared to the genome and vice versa')
parser.add_argument('-i', '--reads', metavar='', required=True, help='the file which contains all the reads in fasta format')
parser.add_argument('-g', '--genome', metavar='', required=True, help='the file that contains the genome')
parser.add_argument('-t', '--threshold', metavar='', required=True, type=int, help='the threshold to filter the kmer')
parser.add_argument('-k', '--size_k', metavar='', required=True, type=int, help='the size of the kmer')

args = parser.parse_args()

file_1 = args.reads  
file_2 = args.genome 
threshold = args.threshold
k = args.size_k

# Function reverse complement which is widely used 

rev = str.maketrans("ACGT", "TGCA")
def reverse_seq(seq: str) -> str:
    return seq.translate(rev)[::-1] 

# ============================================================================ #
#                      First Module : K-mers Dictionary                        #
# ============================================================================ #

def read_file(file):
    """Creation of canonical solid k-mers dictionary

    A dictionary (kmer_dict) having for key the k-mers and for value their abundance 
    is created.
    For each read we create k-mers of size k and their reverse complement. They
    are both compared and only the minimal lexicographical one is retained. If 
    this one not yet in the dictionary it's added, else we adding 1 to its value.  

    Once the dictionary is completed, we remove all the k-mers whose value is 
    lower than t.

    Args:
        file (fasta): file which contains all the read in fasta format
    
    Return:
        dict: canonical solid k-mers dictionary 
    """
    with open(file,"r") as fh: 
        kmer_dict = {}
        count = 0
        for line in fh:
            if ">" not in line:
                for index in range(0,len(line)-k):
                    count += 1 
                    first = line[index:index+k]
                    second = reverse_seq(first)
                    canon = min(first,second)

                    if canon not in kmer_dict:
                        kmer_dict[canon] = 1 
                    else:
                        kmer_dict[canon] += 1

    for kmer, occ in list(kmer_dict.items()):
        if occ < threshold:
            del kmer_dict[kmer]
            
    print(f"Read the {count} {k}-mers from {file}.")
    print(f"Among them {len(kmer_dict)} are solid (abundance >= {threshold}).")
    
    return kmer_dict

# ============================= First Module End ============================= #

# ============================================================================ #
#                      Second Module : SPSS Construction                       #
# ============================================================================ #

def extend_SPSS(kmer_seed, kmer_canon_dict, k2):
    """SPSS's extension creation 

    We create right extension of the interest k-mer using De Bruijn graph. To 
    create the k-mers sons, we remove the first character of the k-mer and we 
    add either A, T, G or C at the end of the same k-mer. We get its canonical
    form and we search it in the dictionary created before. If we found it, we
    remove canonical form to the dictionary and we add the new k-mer to the 
    sequence. If we don't find we pass to the next k-mer son. We repeat the 
    process until the condition no longer true.   

    Args:
        kmer_seed (str): interest k-mer
        kmer_canon_dict (dict): k-mers dictionary 
        k2 : k-mer size

    Returns:
        str: Right extension 
        dict: New k-mer dictionary
    """
    while True:

        first = min(kmer_seed[-k2:] + "A", reverse_seq(kmer_seed[-k2:] + "A"))
        second = min(kmer_seed[-k2:] + "T", reverse_seq(kmer_seed[-k2:] + "T"))
        third = min(kmer_seed[-k2:] + "G", reverse_seq(kmer_seed[-k2:] + "G"))
        fourth = min(kmer_seed[-k2:] + "C", reverse_seq(kmer_seed[-k2:] + "C"))
        if first in kmer_canon_dict:
            del kmer_canon_dict[first]
            kmer_seed += "A"
        elif second  in kmer_canon_dict:
            del kmer_canon_dict[second]
            kmer_seed += "T"
        elif third in kmer_canon_dict:
            del kmer_canon_dict[third]
            kmer_seed += "G"
        elif fourth  in kmer_canon_dict:
            del kmer_canon_dict[fourth]
            kmer_seed += "C"
        else:
            break 
    
    return kmer_seed, kmer_canon_dict


def SPSS(kmer_canon_dict):
    """Create SPSS

    We choose randomly a key in the dictionary to get the first k-mer. We create
    its right extension thanks to the previous function. 
    The left extension of the same k-mer is performed by doing the right extension
    of its reverse complement and then doing the reverse complement of this second
    extension. When their are no more succesors, the symbol '$' is added at the
    end of the sequence and a new k-mer is chosen. We repeat this process until 
    the dictionary is empty.
    The calculation time is also compute.

    Args:
        kmers_cano_dict (dict): Dictionnary contenaing all canonical solid k-mers

    Returns:
        str: spss sequence with "$"
    """
    size_kmer_canon = len(kmer_canon_dict)
    print("Construct spss")
    k2 = k-1
    spss = ""
    concat_numb = 0
    max_prog_bar = len(kmer_canon_dict)
    with Timer() as total_time:
        while len(kmer_canon_dict) != 0:
            update_progress((max_prog_bar - len(kmer_canon_dict))/ max_prog_bar)
            random.seed(1)
            kmer_seed = random.choice(list(kmer_canon_dict.keys()))
            del kmer_canon_dict[kmer_seed]
            
            rev_kmer_seed = reverse_seq(kmer_seed)
            ext_right, kmer_canon_dict = extend_SPSS(kmer_seed, kmer_canon_dict, k2)
            ext_left, kmer_canon_dict = extend_SPSS(rev_kmer_seed, kmer_canon_dict, k2)

            ext_left_rev = reverse_seq(ext_left[k:])
            ext_tot = ext_left_rev + ext_right
            spss += ext_tot + "$"
            concat_numb += 1
        update_progress(1)

        print(f"nb kmers = {size_kmer_canon}")
        print(f"OUT |SPSS(K)| = {len(spss)}")
        print(f"OUT #SPSS(K) = {concat_numb}")
    total_time.print("OUT TIME_SPSS_CONSTRUCTION = {} seconds")
    
    return spss

# ============================ Second Module End ============================= #

# ============================================================================ #
#                       FM INDEX Part (Burrows Wheeler)                        #
# ============================================================================ #

def get_bwt(s: str, sa: List[int]) -> str:
    """Get the Burrow Wheeler Transform of a sequence `s` \
    given its pre-computer suffix array.

    Args:
        s (str): sequence fom which we compute the bwt
        sa (list of int): sa computed from `s`

    Returns:
        (str): the Burrow Wheeler Transform of s
    """  
    bwt = ''
    for pos in sa:
        if pos == 0:
            bwt += '$'
        else:
            bwt += s[pos - 1]

    return bwt

 
def get_r_n(bwt: str) -> List[int]:
    """Get the rank of each character in the given bwt.
       And return the number of occurrences of each character in bwt.

    Args:
        bwt (str): (any) sequence fom which we compute n
    Returns:
        list of int: rank of each character in bwt
        dict : the number of occurrences of each character in bwt
    """
    r = []
    n = {letter: 0 for letter in 'ATGC$'}
    for letter in bwt:
        n[letter] += 1
        r.append(n[letter])

    return r, n


def left_first(alpha: str, k: int, n: Dict[str, int]) -> int:
    """Returns the line of the k^th occurrence of character alpha in the F.

    Args:
        alpha (str): concerned character
        k (int): k^th occurrence of alpha
        n (dict): line of the first occurrence of alpha

    Raises:
        KeyError: if alpha not in alphabet

    Returns:
        int: line corresponding to the k^th occurrence
            of character alpha in F
    """
    if alpha == '$':
        return 0
    if alpha == 'A':
        return n['$'] + k - 1
    if alpha == 'C':
        return n['$'] + n['A'] + k - 1
    if alpha == 'G':
        return n['$'] + n['A'] + n['C'] + k - 1
    if alpha == 'T':
        return n['$'] + n['A'] + n['C'] + n['G'] + k - 1
    raise KeyError


def contains(p: str, bwt: str,
             n: Dict[str, int], r: List[int]) -> Tuple[bool, int, int]:
    """Return is p in original sequence?

    Args:
        p (str): pattern to search
        bwt (str): BWT
        n (dict): line of the first occurrence of alpha
        r (list of int): rank of each character in bwt

    Returns:
        bool: is p in original sequence?
        int: start interval (if True, else -1)
        int: end interval (if True, else -1)
    """
    start = 0
    stop = len(bwt) - 1

    for i in range(len(p) - 1, -1, -1):
        # Find interval beginning
        new_start = first_line_in_bwt_containing_alpha(
            bwt, start, stop, p[i],
        )
        if new_start > stop:  
            return False

        # Find interval ending
        new_stop = last_line_in_bwt_containing_alpha(
            bwt, stop, start, p[i],
        )
        start = left_first(p[i], r[new_start], n)
        stop = left_first(p[i], r[new_stop], n)
        
    return True


def first_line_in_bwt_containing_alpha(bwt: str, start: int, max_line: int,
                                       alpha: str) -> int:
    """Return first BWT line index containing alpha.

    Args:
        bwt (str): BWT
        start (int): start index
        max_line (int): max index
        alpha (str): character to look for

    Returns:
        int: first BWT line index containing alpha, else max_line + 1
    """
    line = start
    while line < max_line + 1 and bwt[line] != alpha:
        line += 1

    return line


def last_line_in_bwt_containing_alpha(bwt: str, stop: int, min_line: int,
                                      alpha: str) -> str:
    """Return last BWT line index containing alpha.

    Args:
        bwt (str): BWT
        stop (int): stop index
        min_line (int): min index
        alpha (str): character to look for

    Returns:
        int: last BWT line index containing alpha, else min_line
    """
    line = stop
    while line > min_line and bwt[line] != alpha:
        line -= 1

    return line

# ============================ FM INDEX Part End ============================= #

# ============================================================================ #
#                       Third Module : FN and FP count                         #
# ============================================================================ #

def FP_index(file_genome, SPSS):
    """Count the number of FP by using FM-index

    We use Burrows-Wheeler Transform to find k-mer from SPSS in the reference 
    genome. If the k-mer is not find we add 1 to the counter. All k-mer with a 
    '$' are automatically dismissed. And this we calculate the rate of FP. 
    The calculation time is also compute.

    Args:
        file_genome (fasta): file which contains the reference genome
        SPSS (str): SPSS sequence
    
    Return:
        int: FP rate             
    """
    with open(file_genome) as fh:
        for line in fh:
            if ">" not in line:
                sa = simple_kark_sort(line+'$')
                bwt = get_bwt(line+'$', sa)
                r,n = get_r_n(bwt)
                count = 0 

                size_seq = len(SPSS)
                iter_max = size_seq - k + 1
                
                with Timer() as total_time: # time all instructions in the ’with’ statements
                    for i in range(iter_max):
                        if i % 100 == 0:
                            update_progress(i/ iter_max)
                        kmer = SPSS[i:i+k]
                        rev_kmer = reverse_seq(kmer)
                        if ('$' not in kmer) and not contains(kmer, bwt, n, r) \
                            and not contains(rev_kmer, bwt, n, r):  
                            count+=1 
                    update_progress(1)
    
    total_time.print("OUT TIME_FP_WITH_INDEX= = {} seconds")
 
    print(f"OUT #FP = {count}")
    fp = (count)*100/iter_max

    return fp 


def FP_in(file_genome, SPSS):
    """Count the number of FP by using FM-index

    We use python pattern matching to find k-mer from SPSS in the reference 
    genome. If the k-mer is not find we add 1 to the counter. All k-mer with a 
    '$' are automatically dismissed. And this we calculate the rate of FP.
    The calculation time is also compute.

    Args:
        file_genome (fasta): file which contains the reference genome
        SPSS (str): SPSS sequence
    
    Return:
        int: FP rate             
    """
    with open(file_genome) as fh:
        for line in fh:
            if ">" not in line:
                count = 0
                index = 0

                size_seq = len(SPSS)
                iter_max = size_seq - k + 1
    
                with Timer() as total_time: 
                    for i in range(iter_max):
                        if i % 100 == 0:
                            update_progress(i/ iter_max)
                        kmer = SPSS[i:i+k]
                        rev_kmer = reverse_seq(kmer)
                        if ('$' not in kmer) and kmer not in line+'$' \
                            and rev_kmer not in line+'$':  
                            count+=1 
                    update_progress(1)
    
    total_time.print("OUT TIME_FP_WITHOUT_INDEX = {} seconds")
    
    print(f"OUT #FP = {count}")
    fp = (count)*100/iter_max

    return fp 


def FN_index(file_genome, SPSS):
    """Count the number of FN by using FM-index

    We use Burrows-Wheeler Transform to find k-mer from the reference genome in
    the SPSS. If the k-mer is not find we add 1 to the counter. All k-mer with 
    a '$' are automatically dismissed. And this we calculate the rate of FN. 
    The calculation time is also compute.

    Args:
        file_genome (fasta): file which contains the reference genome
        SPSS (str): SPSS sequence
    
    Return:
        int: FN rate             
    """
    with open(file_genome) as fh:
        for line in fh:
            if ">" not in line:
                sa = simple_kark_sort(SPSS)
                bwt = get_bwt(SPSS, sa)
                r,n = get_r_n(bwt)
                count = 0 

                size_gen = len(line+'$')
                iter_max = size_gen - k + 1

                with Timer() as total_time: 
                    for i in range(iter_max):
                        if i % 100 == 0:
                            update_progress(i/ iter_max)
                        kmer = line[i:i+k]
                        rev_kmer = reverse_seq(kmer)
                        if (not contains(kmer, bwt, n, r)) \
                            and (not contains(rev_kmer, bwt, n, r)) : 
                            count+=1 
                    update_progress(1)
    
    total_time.print("OUT TIME_FN_WITH_INDEX = {} seconds")
 
    print(f"OUT #FN = {count - 1}") # Because of the last k-mer which contains the last '$' 
    fn = (count - 1)*100/iter_max
    
    return fn 


def FN_in(file_genome, SPSS):
    """Count the number of FN by using FM-index

    We use python pattern matching to find k-mer from reference genome in the 
    SPSS. If the k-mer is not find we add 1 to the counter. All k-mer with a 
    '$' are automatically dismissed. And this we calculate the rate of FN.
    The calculation time is also compute.

    Args:
        file_genome (fasta): file which contains the reference genome
        SPSS (str): SPSS sequence
    
    Return:
        int: FN rate             
    """
    with open(file_genome) as fh:
        for line in fh:
            if ">" not in line:
                count = 0
                index = 0

                size_gen = len(line+'$')
                iter_max = size_gen - k + 1
    
                with Timer() as total_time:
                    for i in range(iter_max):
                        if i % 100 == 0:
                            update_progress(i/ iter_max)
                        kmer = line[i:i+k]
                        rev_kmer = reverse_seq(kmer)
                        if kmer not in SPSS and rev_kmer not in SPSS:  
                            count+=1 
                    update_progress(1)
    total_time.print("OUT TIME_FN_WITHOUT_INDEX = {} seconds")
    
    print(f"OUT #FN = {count - 1}")
    fn = (count - 1)*100/iter_max

    return fn 

# ============================ Third Module End ============================== #

# ============================================================================ #
#                                    Output                                    #
# ============================================================================ #

def main():
    #### Compute Ditionary and SPSS ####
    kmer_dict = read_file(file_1)
    spss = SPSS(kmer_dict)
    print("\n")

    print("Count FP and FN using python pattern matching")
    print("Count kmers from the genome not in the spss")
    error_rate = FN_index(file_2, spss)
    print(f"OUT ERROR_RATE : {round(error_rate, 3)} %")
    print("\n")

    print("Count kmers from the spss not in the genome")
    error_rate = FP_index(file_2, spss)
    print(f"OUT ERROR_RATE : {round(error_rate, 3)} %")
    print("\n")

    print("Count FP and FN using python pattern matching")
    print("Count kmers from the genome not in the spss")
    FN_in(file_2,spss)
    print(f"OUT ERROR_RATE : {round(error_rate, 3)} %")
    print("\n")

    print("Count kmers from the spss not in the genome")
    FP_in(file_2,spss)
    print(f"OUT ERROR_RATE : {round(error_rate, 3)} %")


if __name__ == '__main__':
    main()

# =============================== Output End ================================= #

# ============================== Project End ================================= #