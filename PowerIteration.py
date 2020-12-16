import numpy as np
def diagonal_eigval(eigen_val):
        arr = [[0 for j in range(len(eigen_val))] for i in range(len(eigen_val))]
        x=0
        for i in eigen_val:
            arr[x][x]=pow(i,20)
            x+=1
            if x>3:
                break
        return arr

def Google_Matrix(A, m):
    N = A.shape[0]
    v = np.ones(N)
    
    # Calculate the degree of each node
    KT = np.dot(A.T, v)

    # Normalize the columns
    for i in range(N):
        A.T[i] = A.T[i]/KT[i]
    
    # Add random links
    S = np.ones((N, N))/N
    G = (1-m)*A+m*S

    return G

def Power_Method(G, iter):
    N = G.shape[0]
    x0 = np.ones(N)/N

    for i in range(iter):
        x0 = np.dot(G, x0)

    return x0

n = int(input("Enter the number of rows:"))
m = int(input("Enter the number of columns:"))
print("Enter the entries rowise (separated by space): ")
matrix=[]
for i in range(n):
    a=[]
    for j in range(m):
        a.append(float(input()))
    matrix.append(a)

A_Matrix=np.array(matrix)
# gives adjacency matrix when M is kept 0 in google matrix
# print(Google_Matrix(A_Matrix,0))
print("Result of Power Iteration Algorithm: ")
print(Power_Method(A_Matrix,100))