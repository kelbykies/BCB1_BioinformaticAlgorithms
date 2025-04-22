# Basic Example

#sequence
A = "CGTAATCGTA"
n = len(A)
word = "CG"
w = len(word)

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

# for i in range(0,(pow(4,w) - 1)):
#     Buck[i] = 0

c = 0

for i in range(0,w-1):
    c = c * 4 + D[A[i]]

for j in range(0, n-w+1):
    c = (c*4 + D[A[j+w-1]]) % pow(4,w)
    List[j] = buck[c]
    buck[c] = j + 1

print("Buck", buck)

print("List", List)