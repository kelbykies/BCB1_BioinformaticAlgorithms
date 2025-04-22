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
    if index1 < len(seq1) and index2 < len(seq2):
        if seq1[index1 - 1] == seq2[index2 - 1]:
            # print(seq1[index1 - 1], seq2[index2 - 1])
            if seq1[index1] == seq2[index2]:
                return matchScore + bonus
            return matchScore
        else:
            return mismatchScore
    else:
        if seq1[index1 - 1] == seq2[index2 - 1]:
            return matchScore
        else:
            return mismatchScore

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

# seq1 = "GATCGTAGAGTGAGACCTAGTGTTTG"
# seq2 = "CTCGTAGGTGAGATTCCTAGTGCC"
seq1 = "ATG"
seq2 = "ACG"
# seq1 = "ATCGGGGGGGGGGG"
# seq2 = "TTTTTTTTTTTATC"


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
    if mat == "S":
        if i == m or j == n or S_matrix[i][j] == 0:

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
# print(row_first, col_first)
# print(row_last, col_last)
#
for i in range(len(OA)):
    if OA[i] == OB[i]:
        match = match + 1
    elif OA[i] == "-" or OB[i] == "-":
        gap_length = gap_length + 1
    else:
        mismatch = mismatch + 1

align_length = len(OA)
percent_id = (int)(100 * (match/align_length))

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
# # print("Length of Alignment: ", align_length)
# # print("Percent identity: ", percent_id, "%")
# print("Number of Matches:", match_count)
# print("Number of Mismatches:", mismatch_count)
# print("Total Length of Gaps:", gap_length)
print("Start Position in A:", rowfirst)
print("Start Position in B:", colfirst)
print("End Position in A:", rowlast - 1)
print("End position in B:", collast - 1, "\n")
#
#
# # Alignment in sections of 70 chars, with each section consisting of 3 rows
print(OA)
print(middle)
print(OB)
