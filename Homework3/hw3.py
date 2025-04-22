from sys import argv


def make_reverse_complement(seq):
    """ Takes any give sequence and returns the reverse complement. """

    reverse_seq = seq[::-1]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_complement = "".join(complement.get(base, base) for base in reverse_seq)

    return rev_complement


def extract_info(input_file):
    """ Reads in the read sequences from the provided sequence file. Sequences and their Reverse Compliments are
    concatenated together.
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
# script, seq_file, wm_file, wlcut = argv
#
# wm_file = open(wm_file, "r")
# word_model = wm_file.read()

# Length of superword
# wlcut = int(wlcut)
# sequence = extract_info(seq_file)

sequence = "ACGT#ACGT"
word_model = "1"
wlcut = 1
checked = word_model.count("1")
unchecked = word_model.count("0")

if checked == 0:
    print("Error. Word Model must contain at least 1 checked position.")
    quit()

w = checked + unchecked
n = len(sequence)


# Computing the code array
index = 0
c = 0
word_code = [0] * n

# Dictionary with base values
D = {}
D['A'] = 0
D['C'] = 1
D['G'] = 2
D['T'] = 3

while index < n - w + 1:
    wm_index = 0
    c = 0
    for j in range(index, index + w):
        if j > n - 1:
            c = -1
            break

        if sequence[j] not in {'A', 'T', 'C', 'G'}:
            c = -1
            break

        else:
            if word_model[wm_index] == "1":
                # if sequence[j] not in {'A', 'T', 'C', 'G'}:
                #     c = -1
                #     break
                c = (c*4 + D[sequence[j]]) % pow(4,w)
            wm_index = wm_index + 1

    word_code[index] = c
    index = index + 1

for i in range(n - w + 1, n):
    word_code[i] = -1

print(word_code)
# Computation of Superword Array
SW = [0] * n
List = [0] * (n + 2)

for i in range(0, n):
    SW[i] = i + 1

buck_size = 1
for i in range(1, checked + 1):
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

    for k in range(j + 1, n + 1):
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