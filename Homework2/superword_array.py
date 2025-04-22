seq = "CATCGTGANATCGTGTN"
n = len(seq)
w = 2
wlcut = 3

#Look up table
buck = [0] * pow(4,w)
List = [0] * n
# Buck = []
# Hash table
D = {}
D['A'] = 0
D['C'] = 1
D['G'] = 2
D['T'] = 3

# Length of word
w = 2

# Length of superword
wlcut = 3

#Look up table
buck = [0] * pow(4,w)
List = [0] * n
word_code = [0] * n


# Hash table
D = {}
D['A'] = 0
D['C'] = 1
D['G'] = 2
D['T'] = 3
c = 0
last = 0

for i in range(0, w):
    if seq[i] not in {'A','T','C','G'}:
        last = i
        if i - last < w:
            c = -1
    else:
        c = c*4 + D[seq[i]]
word_code[0] = c

for j in range(1, n - w + 1):
    if seq[j + 1] not in {'A','T','C','G'}:
        last = j + 1
    # if seq[j] not in {'A','T','C','G'}:
    #     last = j
    if (j - last) < (w - 1):
        word_code[j] = -1
        c = (c*4 + 0) % pow(4,w)
    else:
        c = (c * 4 + D[seq[j + w - 1]]) % pow(4, w)
        # List[j] = buck[c]
        # buck[c] = j
        word_code[j] = c

word_code[n-1] = -1

sw = [0] * n
sw[0] = 0
# sw[n+1] = 0


for i in range(0,n):
    sw[i] = i + 1

bsize = 1
for i in range(1, w + 1):
    bsize = bsize * 4

# for wlev in range(1, wlcut):
for wlev in range(0, wlcut):
    for c in range(0, bsize):
        buck[c] = 0
    for i in range(1, n + 1):
        # p = sw[i] - w
        p = sw[i - 1] - w
        if p > 0:
            c = word_code[p]
            List[p] = buck[c]
            buck[c] = p
    for p in range(n - w + 1, n + 1):
        # c = word_code[p]
        c = word_code[p - 1]
        # print(c)
        # List[p] = buck[c]
        # buck[c] = p
        List[p-1] = buck[c]
        buck[c] = p
    # k = n + 1
    k = n
    for c in range(bsize - 1, -2, -1):
        p = buck[c]
        while p > 0:
            sw[k - 1] = p
            k = k - 1
            p = List[p - 1]
    buck[c] = k
    # buck[bsize] = n + 1
    buck[bsize - 1] = n + 1

print(sw)