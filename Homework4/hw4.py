# Filename: hw4.py
# Author: Kelby Kies
# Date Created: 11/21/2019
# Date Last Modified: 12/4/2019
# Python Version: Python 3.6

from sys import argv


def make_reverse_complement(seq):
    """ Takes any give sequence and returns the reverse complement. """

    reverse_seq = seq[::-1]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_complement = "".join(complement.get(base, base) for base in reverse_seq)

    return rev_complement


def extract_seq(file, number_of_file):
    """ Extracts the name, sequence and length from a FASTA seq_file. """
    file = open(file)
    file_contents = file.read()
    seq = ""

    for line in file_contents.splitlines():
        if line.startswith(">"):
            key = line.strip()[1:]
        else:
            if number_of_file == 0:
                # Do not include the Reverse Complement on the first file
                seq = seq + line + "#"
            else:
                seq = seq + line + "#" + make_reverse_complement(line) + "#"

    section_length = len(seq)
    file.close()
    return seq, section_length


def create_alignment(genome_array, seq, num):
    """ Composes the final multiple sequence alignment for each genome."""

    file_alignment = []

    for file in range(0,num):
        file_alignment.append("")

    for a in range(0, num):
        current = []

        for b in range(0, len(genome_array)):

            tup = (genome_array[b][a], seq[genome_array[b][a]:genome_array[b][a] + (wlcut * w)])

            current.append(tup)

            current.sort(key=lambda x: x[0])

        for d in current:

            file_alignment[a] = file_alignment[a] + d[1]

    return file_alignment


# Reading in the paramters
script, file_of_filenames, wm_file, wlcut = argv

wm_file = open(wm_file, "r")
word_model = wm_file.read()

# Length of superword
wlcut = int(wlcut)

# List of Files
read_files = open(file_of_filenames, "r")
file_names = read_files.readlines()
m = len(file_names)

# Compose full sequence by concatenating sequences from each file
st = [0] * m
section_start = 0
final_seq = ""
for i in range(0, m):
    new_file = file_names[i]
    st[i] = section_start
    seq_part = extract_seq(new_file[:-1], i)
    final_seq = final_seq + seq_part[0]
    section_start = section_start + seq_part[1]


checked = word_model.count("1")
unchecked = word_model.count("0")
n = len(final_seq)
w = checked + unchecked

if checked == 0:
    print("Error. Word Model must contain at least 1 checked position.")
    quit()

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

        if final_seq[j] not in {'A', 'T', 'C', 'G'}:
            c = -1
            break

        else:
            if word_model[wm_index] == "1":
                if final_seq[j] not in {'A', 'T', 'C', 'G'}:
                    c = -1
                    break
                c = (c*4 + D[final_seq[j]]) % pow(4,w)
            wm_index = wm_index + 1

    word_code[index] = c
    index = index + 1

for i in range(n - w + 1, n):
    word_code[i] = -1

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

overlap_array = list('n' * n)
mult_seq = ""
block_array = [0] * m
j = 0
num_of_blocks = 0
genome_array = []

while j < n:
    s = []
    x = 1
    for wlev in range(0, wlcut):
        one = SW[j] + (wlev * w)
        if one >= n:
            x = 0
            break
        else:
            if word_code[one - 1] == -1:
                x = 0
                break

    if x == 0:
        j = j + 1
        continue

    superword = final_seq[SW[j] - 1: SW[j] - 1 + (wlcut * w)]
    s.append(SW[j] - 1)

    for k in range(j + 1, n + 1):
        y = 1
        if k == n:
            break
        else:
            for wlev in range(0, wlcut):
                one = SW[j] + (wlev * w)
                two = SW[k] + (wlev * w)
                if two >= n:
                    y = 0
                    break
                if word_code[two - 1] == -1:
                    y = 0
                    break

                if word_code[one - 1] != word_code[two - 1]:
                    y = 0
                    break
                else:
                    s.append(SW[k] - 1)

            if y == 0:
                break

    block_length = k - j
    # Check if Block Length is equal to m
    if block_length == m:
        block_position = j
        num = 0

        # Check to see if SW comes from Unique Genomes
        for y in range(0, m):

            if y < m - 1:
                genome = final_seq[st[y]:st[y + 1]]
                if genome.count(superword) == 1:
                    num = num + 1

            else:
                genome = final_seq[st[y]: n]
                if genome.count(superword) == 1:
                    num = num + 1

            block_position = block_position + 1

        if num == m:
            # Check for Overlap
            count = 0
            for p in range(j, k):
                overlap = overlap_array[SW[p] - 1: SW[p] - 1 + (wlcut * w)]

                if "y" not in overlap:
                    count = count + 1
                    overlap_array[SW[p] - 1: SW[p] - 1 + (wlcut * w)] = 'y' * (wlcut * w)

            if count == m:
                num_of_blocks = num_of_blocks + 1
                genome_array.append(s)

    j = k


subset_length = num_of_blocks * (w * wlcut)
mult_seq = create_alignment(genome_array, final_seq, m)


# Print the alignment summary
summary_name = file_of_filenames + ".summary"
fh = open(summary_name, "w")
fh.write("Sequence of Word Model: ")
fh.write(word_model)
fh.write("\n")
fh.write("Value for wlcut: ")
fh.write(str(wlcut))
fh.write("\n")
fh.write("Length of largest subset of SW blocks: ")
fh.write(str(subset_length))
fh.write("\n")
fh.write("# of SW blocks in Subset: ")
fh.write(str(int(num_of_blocks)))
fh.close()

# Print Alignment to a file in FASTA format
align_name = file_of_filenames + ".alignment"
f = open(align_name, "w")

for i in range(0, m):
    if len(mult_seq) > 0:
        f.write(">")
        f.write(file_names[i])
        f.write(mult_seq[i] + "\n")
    else:
        f.write(">")
        f.write(file_names[i])
        f.write("\n")
f.close




