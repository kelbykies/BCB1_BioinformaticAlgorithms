# Filename: hw1.py
# Author: Kelby Kies
# Date Created: 9/9/2019
# Date Last Modified: 9/23/2019
# Python Version: Python 3.6


# Imports
from sys import argv


def extract_seq(file):
    """ Extracts the name, sequence and length from a FASTA seq_file. """
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
    return name, seq, length


def make_matrix(column, row):
    """ Creates a 2D matrix """
    matrix = [0] * column
    for z in range(column):
        matrix[z] = [0] * row
    return matrix


def determine_match(index1, index2, seq1, seq2, bonus):
    """
    Determines if a nucleotide in the first sequence matches/mismatches a nucleotide in the second sequence.
    If the match_count is consecutive then the bonus is added to the match_count score.
    """
    if index1 < len(seq1) and index2 < len(seq2):
        if seq1[index1 - 1] == seq2[index2 - 1]:
            if seq1[index1] == seq2[index2]:
                return match_score + bonus
            return match_score
        else:
            return mismatch_score
    else:
        if seq1[index1 - 1] == seq2[index2 - 1]:
            return match_score
        else:
            return mismatch_score


def print_alignment(align1, middle, align2, index1, index2):
    """
    Prints sequence alignment in sections of 70 chars. Each section consists of 3 rows.
    Each row begins with the position of the first nucleotide in the sequence.

    """
    start = 0
    end = start + 70
    output = ""
    offset = 10
    while start <= len(align1):
        substring = str(index1) + " " * (offset - len(str(index1))) + align1[start:end] + "\n" \
                    + " " * offset + middle[start:end] + "\n"\
                    + str(index2) + " " * (offset - len(str(index2))) + align2[start:end] + "\n"
        output = output + substring + "\n"
        start = end
        end = start + 70
        index1 = index1 + 70
        index2 = index2 + 70
    print(output)


# Alignment scoring parameters
match_score = 10
mismatch_score = -15
gap_open = 40
gap_extension = 2

score = 0
bonus = 0
match_count = 0
mismatch_count = 0
align_length = 0
percent_id = 0
gap_length = 0
m = 0
n = 0

# Input from user
script, file1, file2, bonus = argv

bonus_score = int(bonus)

seq1_name = extract_seq(file1)[0]
seq2_name = extract_seq(file2)[0]
seq1 = extract_seq(file1)[1]
seq2 = extract_seq(file2)[1]

m = extract_seq(file1)[2]        # rows
n = extract_seq(file2)[2]        # columns

row_first = m
col_first = n

# Create Matrices: S, D, I
s_matrix = make_matrix(m + 1, n + 1)
d_matrix = make_matrix(m + 1, n + 1)
i_matrix = make_matrix(m + 1, n + 1)

# Compute Matrices
s_matrix[m][n] = 0
d_matrix[m][n] = -(gap_open + gap_extension)
i_matrix[m][n] = -(gap_open + gap_extension)

for j in range(n - 1, -1, -1):
    s_matrix[m][j] = 0
    i_matrix[m][j] = -(gap_open + gap_extension)
    d_matrix[m][j] = -(gap_open + gap_extension)
for i in range(m - 1, -1, -1):
    # outer loop
    s_matrix[i][n] = 0
    d_matrix[i][n] = -(gap_open + gap_extension)
    i_matrix[i][n] = -(gap_open + gap_extension)
    for j in range(n - 1, -1, -1):
        # inner loop
        d_matrix[i][j] = max(d_matrix[i + 1][j] - gap_extension, s_matrix[i + 1][j] - gap_open - gap_extension)
        i_matrix[i][j] = max(i_matrix[i][j + 1] - gap_extension, s_matrix[i][j + 1] - gap_open - gap_extension)
        s_matrix[i][j] = max(0, s_matrix[i + 1][j + 1] + determine_match(i + 1, j + 1, seq1, seq2, bonus_score),
                             d_matrix[i][j], i_matrix[i][j])

        if score < s_matrix[i][j]:
            score = s_matrix[i][j]
            row_first = i
            col_first = j


# Traceback
opt_align1 = ""
middle = ""
opt_align2 = ""
i = row_first
j = col_first
matrix = "S"

while i <= m and j <= n:
    if matrix == "S":
        if i == m or j == n or s_matrix[i][j] == 0:
            break

        if s_matrix[i][j] == d_matrix[i][j]:
            matrix = "D"
            continue

        if s_matrix[i][j] == i_matrix[i][j]:
            matrix = "I"
            continue
        opt_align1 = opt_align1 + seq1[i]
        middle = middle + "|"
        opt_align2 = opt_align2 + seq2[j]
        i = i + 1
        j = j + 1
        continue

    if matrix == "D":
        opt_align1 = opt_align1 + seq1[i]
        middle = middle + "-"
        opt_align2 = opt_align2 + " "
        if (i == m - 1) or (d_matrix[i][j] == s_matrix[i + 1][j] - gap_open - gap_extension):
            matrix = "S"
        i = i + 1
        continue

    if matrix == "I":
        opt_align1 = opt_align1 + " "
        middle = middle + "-"
        opt_align2 = opt_align2 + seq2[j]
        if (j == n - 1) or (i_matrix[i][j] == s_matrix[i][j + 1] - gap_open - gap_extension):
            matrix = "S"
        j = j + 1
        continue

row_last = i + 1
col_last = j + 1
row_first = row_first + 1
col_first = col_first + 1

for x in range(len(opt_align1)):
    if opt_align1[x] == opt_align2[x]:
        match_count = match_count + 1
    elif opt_align1[x] != opt_align2[x] and middle[x] != "-":
        mismatch_count = mismatch_count + 1
    else:
        gap_length = gap_length + 1

align_length = len(opt_align1)
percent_id = int(100 * (match_count / align_length))

# Summary
print("Scoring Parameters:\n")
print("Match Score:", match_score)
print("Mismatch Score:", mismatch_score)
print("Gap-Open Penalty:", gap_open)
print("Gap-Extension Penalty:", gap_extension)
print("Bonus Parameter: ", bonus,  "\n")

print("Sequence 1 Name:", seq1_name)
print("   Length:", m, "\n")
print("Sequence 2 Name:", seq2_name)
print("   Length:", n, "\n")

print("Alignment Score:", score)
print("Length of Alignment: ", align_length)
print("Percent identity: ", percent_id, "%")
print("Number of Matches:", match_count)
print("Number of Mismatches:", mismatch_count)
print("Total Length of Gaps:", gap_length)
print("Start Position in A:", row_first)
print("Start Position in B:", col_first)
print("End Position in A:", row_last)
print("End position in B:", col_last, "\n")

print_alignment(opt_align1, middle, opt_align2, row_first, col_first)



