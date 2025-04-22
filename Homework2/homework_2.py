from sys import argv


def make_reverse_complement(seq):
    """

    """

    reverse_seq = seq[::-1]
    rev_complement = ""

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_complement = "".join(complement.get(base, base) for base in reverse_seq)

    return rev_complement


def extract_info(input_file):
    """


    """

    seq = ""
    file = open(input_file, "r")

    line_number = 1
    seq_number = 2
    for line in file:
        if line_number % seq_number == 0:

            seq = seq + line[:-1] + "#" + make_reverse_complement(line[:-1]) + "#"
            seq_number = seq_number + 4
        line_number = line_number + 1

    file.close()
    return seq


# Reading in the FASTQ File
# script, seq_file, wlcut = argv
# seq_file = seq_file
# #Length of superword
# wlcut = int(wlcut)
# seq = extract_info(seq_file)


# output = open("output2.txt", "w+")
# seq = "CGTAATCGTA#"
seq = "AGAG#"
wlcut = 2
# # Length of sequence
# seq = "TT#"
n = len(seq)

# Length of word
w = 1

word_code = [0] * n

# Computation of Code
D = {}
D['A'] = 0
D['C'] = 1
D['G'] = 2
D['T'] = 3

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
            c = (c*4 + D[seq[j]]) % pow(4, w)

    word_code[i] = c
    i = i + 1
word_code[n - 1] = -1
print(word_code)

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

# for i in range(len(SW)):
#
#     output.write(str(SW[i])+"\n")

# Computing the frequency
freq = [0] * (n + 1)
j = 0
while j < n:
    x = 1

    for wlev in range(0, wlcut):
        one = SW[j] + wlev * w
        if one > n:
            x = 0
            break
        else:
            if word_code[one - 1] == -1:
                x = 0
                break

    if x == 0:
        j = j + 1
        continue

    for k in range(j+1, n + 1):
        y = 1
        if k == n:
            break
        else:
            for wlev in range(0, wlcut):
                one = SW[j] + wlev * w
                two = SW[k] + wlev * w

                if two >= n:
                    y = 0
                    break
                if word_code[two - 1] == -1:
                    y = 0
                    break

                if word_code[one - 1] != word_code[two - 1]:
                    y = 0
                    break

            if y == 0:
                break

    freq[k - j - 1] = freq[k - j - 1] + 1
    j = k

print("k: count(k):")
for i in range(0, n):
    if freq[i] != 0:
        print(i + 1, "  ", freq[i])

