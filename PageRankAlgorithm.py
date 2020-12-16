''' NOTE: The Pagerank Algorithm is essentially not a part of our defination. It has been used to illustrate the fact steady state is actually,
# what we intended to find through Page Rank Algorithm as well. The page rank algorithm has been imported from Wikepedia.'''

import numpy as np

def pagerank(M, num_iterations: int = 4500, d: float = 0.85):
    """PageRank: The trillion dollar algorithm."""
    N = M.shape[1]
    v = np.random.rand(N, 1)
    v = v / np.linalg.norm(v, 1)
    M_hat = (d * M + (1 - d) / N)
    for i in range(num_iterations):
        v = M_hat @ v
    return v

              
n = int(input("Enter the number of rows:"))
m = int(input("Enter the number of columns:"))
print("Enter the entries rowise (separated by space): ")
matrix=[]
for i in range(n):
    a=[]
    for j in range(m):
        a.append(float(input()))
    matrix.append(a)

M=np.array(matrix)

# v = pagerank(M, 100, 0.86)
v_absolute=pagerank(M, 100, 1)
print("Pagerank of given Network: ")
print(v_absolute)
# print("Real time Value: ")
# print(v)