from sys import argv
from Bio import SeqIO

# def extract(seq_file):
#     for seq_record in SeqIO.parse(seq_file, "fasta"):
#         name = seq_record.id
#         seq = seq_record.seq
#         length = len(seq_record)
#     return name, seq, length
#
#
# script, file1, file2, bonus = argv
# for seq_record in SeqIO.parse(file1, "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     # print(len(seq_record))
#
# seq1_name = extract(file1)[0]
# seq2_name = extract(file2)[0]
# seq1 = extract(file1)[1]
# seq2 = extract(file2)[1]
# m = extract(file1)[2]
# n = extract(file2)[2]

# print(n)
#
def extract_seq(file):

    file = open(file)
    file_contents = file.read()
    seq = ""
    for line in file_contents.splitlines():
        if line.startswith(">"):
            key = line.strip()[1:]
            name = key.split()[0]
        else:
            seq = seq + line
    length = len(seq)
    file.close()
    # print(length)
    return name, seq, length


script, file1, file2, bonus = argv
bonus_score = bonus
extract_seq(file2)
seq1_name = extract_seq(file1)[0]
seq2_name = extract_seq(file2)[0]
seq1 = extract_seq(file1)[1]
seq2 = extract_seq(file2)[1]
m = extract_seq(file1)[2]
n = extract_seq(file2)[2]

print(seq2_name, n)
