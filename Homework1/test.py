import sys
from sys import argv
import textwrap


def make_matrix(column, row):
    """
    :param m:
    :param n:
    :return:
    """
    matrix = [-1] * column
    for i in range(column):
        matrix[i] = [-1] * row
    return matrix


def determine_match(index1, index2, seq1, seq2, bonus):
    """
    :param seq_1_index:
    :param seq_2_letter_index:
    :param seq1:
    :param seq2:
    :return:
    """
    # print(index1, index2, len(seq1), len(seq2))
    # if index1 < len(seq1) and index2 < len(seq2):
    #     if seq1[index1] == seq2[index2]:
    #         if seq1[index1 - 1] == seq2[index2 - 1]:
    #             return match_score + bonus
    #         else:
    #             return match_score
    #     else:
    #         return mismatch_score
    # else:
    if seq1[index1 - 1] == seq2[index2 - 1]:
        return matchScore
    else:
        return mismatchScore


def print_alignment(a1, a2, a3, index1, index2):
    start = 0
    end = start + 70
    output = ""
    offset = 10
    while start <= len(a1):
        # substring = str(index1) + " " * offset + a1[start:end] + "\n" + " " * offset + a2[start:end]\
        #             + "\n" + str(index2) + " " * offset + a3[start:end] + "\n"
        substring =  str(index1) + " " * (offset - len(str(index1))) + a1[start:end] + "\n" + " " * offset + a2[start:end] + "\n" + str(index2) + " " * (offset - len(str(index1))) + a3[start:end] + "\n"
        output = output + substring + "\n"
        start = end
        end = start + 70
        index1 = index1 + 70
        index2 = index2 + 70
    return output


# Alignment scoring parameters
matchScore = 10
# mismatch_score = -15
mismatchScore = -20
gapOpen = 40        #q
gapExtension = 2    #r
score = 0
bonus = 1
match = 0
mismatch = 0
align_length = 0
percent_id = 0
gap_length = 0

seq1 = "GTCGTAGAGTGAGACCTAGTGTTTG"
seq2 = "CTCGTAGGTGAGATTCCTAGTGCC"
# seq1 = "AG"
# seq2 = "AG"

m = len(seq1)        # m; rows
n = len(seq2)        # n; columns

rowfirst = m
colfirst = n

S_matrix = make_matrix(m + 1, n + 1)
D_matrix = make_matrix(m + 1, n + 1)
I_matrix = make_matrix(m + 1, n + 1)

S_matrix[m][n] = 0
D_matrix[m][n] = -(gapOpen + gapExtension)
I_matrix[m][n] = -(gapOpen + gapExtension)

for j in range(n - 1, -1, -1):
    S_matrix[m][j] = 0
    I_matrix[m][j] = -(gapOpen + gapExtension)
    D_matrix[m][j] = -(gapOpen + gapExtension)
for i in range(m - 1, -1, -1):
    # outer loop
    S_matrix[i][n] = 0
    D_matrix[i][n] = -(gapOpen + gapExtension)
    I_matrix[i][n] = -(gapOpen + gapExtension)
    for j in range(n - 1, -1, -1):
        # inner loop
        D_matrix[i][j] = max(D_matrix[i + 1][j] - gapExtension, S_matrix[i + 1][j] - gapOpen - gapExtension)
        I_matrix[i][j] = max(I_matrix[i][j + 1] - gapExtension, S_matrix[i][j + 1] - gapOpen - gapExtension)
        # s_matrix[i][j] = max(0, s_matrix[i + 1][j + 1] + determine_match(i + 1, j + 1, seq1, seq2, bonus),
        #                      d_matrix[i][j], i_matrix[i][j])
        S_matrix[i][j] = max(0, S_matrix[i + 1][j + 1] + determine_match(i + 1, j + 1, seq1, seq2, bonus),
                             D_matrix[i][j], I_matrix[i][j])

        if score < S_matrix[i][j]:
            score = S_matrix[i][j]
            rowfirst = i
            colfirst = j

# for z in d_matrix:
#     print(z)


# Traceback
OA = ""
middle = ""
OB = ""
# i = row_first
# j = col_first
i = rowfirst
j = colfirst
mat = "S"
match = 0
mismatch = 0
align_length = 0
percent_id = 0
gap_length = 0



while i <= m and j <= n:
    # print("tb:i={}\tj={}\tmatrix={}".format(i, j, matrix))
    # print(opt_align1)
    # print(middle)
    # print(opt_align2)
    if mat == "S":
        if i == m or j == n or S_matrix[i][j] == 0:
            # print(s_matrix[i][j])
            break
        if S_matrix[i][j] == D_matrix[i][j]:
            mat = "D"
            continue
        if S_matrix[i][j] == I_matrix[i][j]:
            mat = "I"
            continue
        OA = OA + seq1[i]
        middle = middle + "|"
        OB = OB + seq2[j]
        # opt_align1 = opt_align1 + seq1[i + 1]
        # opt_align2 = opt_align2 + seq2[j + 1]
        # print(seq1[i-1], seq2[j - 1])
        i = i + 1
        j = j + 1

        continue
    if mat == "D":
        OA = OA + seq1[i]
        # opt_align1 = opt_align1 + seq1[i + 1]
        middle = middle + "-"
        OB = OB + "-"

        if (i == m - 1) or (D_matrix[i][j] == S_matrix[i+1][j] - gapOpen - gapExtension):
            mat = "S"
        i = i + 1

        continue
    if mat == "I":
        OA = OA + "-"
        middle = middle + "-"
        OB = OB + seq2[j]
        # opt_align2 = opt_align2 + seq2[j + 1]
        if (j == n - 1) or (I_matrix[i][j] == S_matrix[i][j + 1] - gapOpen - gapExtension):
            mat = "S"
        j = j + 1

        continue

rowlast = i
collast = j
print(rowfirst, colfirst)
print(rowlast, collast)

for i in range(len(OA)):
    if OA[i] == OB[i]:
        match = match + 1
    elif OA[i] == "-" or OB[i] == "-":
        gap_length = gap_length + 1
    else:
        mismatch = mismatch + 1

# align_length = len(opt_align1)
# percent_id = (int)(100 * (match_count/align_length))

# Summary
print("Scoring Parameters:\n")
print("Match Score:", matchScore)
print("Mismatch Score:", mismatchScore)
print("Gap-Open Penalty:", gapOpen)
print("Gap-Extension Penalty:", gapExtension, "\n")
print("Bonus Parameter: ", bonus)

print("Sequence A:", seq1) #Name
print("Length:", m, "\n")   #Length
print("Sequence B:", seq2)  #Name
print("Length:", n, "\n") # Length

print("Alignment Score:", score)
# print("Length of Alignment: ", align_length)
# print("Percent identity: ", percent_id, "%")
print("Number of Matches:", match)
print("Number of Mismatches:", mismatch)
print("Total Length of Gaps:", gap_length)
print("Start Position in A:", rowfirst)
print("Start Position in B:", colfirst)
print("End Position in A:", rowlast - 1)
print("End position in B:", collast - 1, "\n")


# Alignment in sections of 70 chars, with each section consisting of 3 rows
print(OA)
print(middle)
print(OB)

# print(print_alignment(opt_align1, middle, opt_align2, row_first, col_first))
