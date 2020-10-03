import random
import sys

n = int(sys.argv[1])
p = float(sys.argv[2])
labelmod = int(sys.argv[3])

edges = []
for i in range(n - 1):
    for j in range(i + 1, n):
        if random.random() < p:
            edges.append((i, j))

print(n, len(edges))
print(" ".join(str(i % labelmod) for i in range(n)))

for u, v in edges:
    print(u, v)
