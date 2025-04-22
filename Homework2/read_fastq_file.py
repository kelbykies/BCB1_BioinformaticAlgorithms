from sys import argv
import re
import string

def reverse_complement(seq):
    reverse_seq = seq[::-1]
    rev_complement = ""

    # rev_complement = reverse_seq.translate(trans)
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for base in reverse_seq:
        if base == 'A' or base =='C' or base == 'T' or base =='G':
            rev_complement = rev_complement + ''.join(complement[base])
        else:
            rev_complement = rev_complement + base
    return rev_complement


def extract_info(file):

    file = open(file)
    file_contents = file.read()
    seq = ""
    file_array = file_contents.splitlines()

    for line in file_array:
        if line.startswith("@"):
            # print(line)
            index_of_line = file_array.index(line)
        # #     print(line)
        # if re.match("^[A|C|T|G]", line[0]):
        #     seq = seq + file_array[index_of_line + 1]
            seq = seq + file_array[index_of_line + 1] + "#" + reverse_complement(file_array[index_of_line + 1]) + "#"
    file.close()
    return seq


script, file, wlcut = argv
file = file
wlcut = wlcut

extract_info(file)

