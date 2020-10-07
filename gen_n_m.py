import random
import sys

n = int(sys.argv[1])
m = int(sys.argv[2])
labelmod = int(sys.argv[3])

edges = random.sample([(i, j) for i in range(n-1) for j in range(i+1, n)], m)

print(n, len(edges))
print(" ".join(str(i % labelmod) for i in range(n)))

for u, v in edges:
    print(u, v)
