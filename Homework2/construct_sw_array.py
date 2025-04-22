def reverse_complement(seq):
    reverse_seq = seq[::-1]
    rev_complement = ""

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for base in reverse_seq:
        rev_complement = rev_complement + ''.join(complement[base])

    return rev_complement


def extract_info(file):

    file_contents = file.read()
    seq = ""
    file_array = file_contents.splitlines()

    for i in range(0, len(file_array) - 1, 4):
        if file_array[i].startswith("@"):
        # This is problematic because some lines start with @ but aren't the beginning
            seq = seq + file_array[i + 1] + "#" + reverse_complement(file_array[i + 1])
    return seq

seq = "CGTAATCGTA#"
# seq = "A#"
# seq = "AACACG#TTT##"
# seq = "AG#"
# # Length of sequence
n = len(seq)
# Length of word
w = 1
wlcut = 1
word_code = [0] * n

# Computation of Code
D = {}
D['A'] = 0
D['C'] = 1
D['G'] = 2
D['T'] = 3

c = 0
#Turned in
# for i in range(0, w):
#     if seq[i] not in {'A', 'T', 'C', 'G'}:
#         c = -1
#     c = c*4 + D[seq[i]]
# word_code[0] = c
# for j in range(1, n - w + 1):
#     if seq[j] not in {'A', 'T', 'C', 'G'}:
#         c = -1
#     else:
#         c = (c*4 + D[seq[j]]) % pow(4, w)
#     word_code[j-1] = c
#
# word_code[n-1] = -1

# #Correct
i = 0
while i < n - w + 1:
    c = 0
    for j in range(i, i + w):
        if j > n-1:
            c = -1
            break
        else:
            if seq[j] not in {'A', 'T', 'C', 'G'}:
                c = -1
                break

    if c != -1:
        c = c*4 + D[seq[i]]
        c = (c*4 + D[seq[i+w-1]]) % pow(4,w)
    word_code[i] = c
    i = i + 1
word_code[n - 1] = -1


# Computation of Superword Array
SW = [0] * n
List = [0] * (n + 2)

for i in range(0, n):
    SW[i] = i + 1

buck_size = 1
for i in range(1, w + 1):
    buck_size = buck_size * 4

buck = [0] * (buck_size + 2)
for wlev in range(1, wlcut + 1):

    for c in range(-1, buck_size):
        buck[c] = 0

    for i in range(1, n + 1):
        p = SW[i - 1] - w

        if p > 0:
            c = word_code[p - 1]
            List[p - 1] = buck[c]
            buck[c] = p

    for p in range(n - w + 1, n + 1):
        c = word_code[p - 1]
        List[p - 1] = buck[c]
        buck[c] = p

    # Phase 2
    k = n
    for c in range(buck_size - 1, -2, -1):
        p = buck[c]
        while p > 0:
            SW[k - 1] = p

            p = List[p - 1]
            k = k - 1
        buck[c] = k

    buck[buck_size - 1] = n + 1