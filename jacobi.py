from numpy import array,identity,diagonal
from math import sqrt
import numpy as np



def Jacobi(A,tol = 1.0e-9): # Jacobi method
# Find largest off-diagonal element a[k,l]
 def maxElem(A):
    n = len(A)
    Amax = 0.0
    for i in range(n-1):
        for j in range(i+1,n):
            if abs(A[i,j]) >= Amax:
                Amax = abs(A[i,j])
                k = i; l = j
    return Amax,k,l

# Rotate to make A[k,l] = 0 and define the rotation matrix
 def rotate(A,p,k,l):
    n = len(A)
    Adiff = A[l,l] - A[k,k]
    if abs(A[k,l]) < abs(Adiff)*1.0e-36: t = A[k,l]/Adiff
    else:
        phi = Adiff/(2.0*A[k,l])
        t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
        if phi < 0.0: t = -t
    c = 1.0/sqrt(t**2 + 1.0); s = t*c
    tau = s/(1.0 + c)
    temp = A[k,l]
    A[k,l] = 0.0
    A[k,k] = A[k,k] - t*temp
    A[l,l] = A[l,l] + t*temp
    for i in range(k): # Case of i < k
        temp = A[i,k]
        A[i,k] = temp - s*(A[i,l] + tau*temp)
        A[i,l] = A[i,l] + s*(temp - tau*A[i,l])
    for i in range(k+1,l): # Case of k < i < l
        temp = A[k,i]
        A[k,i] = temp - s*(A[i,l] + tau*A[k,i])
        A[i,l] = A[i,l] + s*(temp - tau*A[i,l])
    for i in range(l+1,n): # Case of i > l
        temp = A[k,i]
        A[k,i] = temp - s*(A[l,i] + tau*temp)
        A[l,i] = A[l,i] + s*(temp - tau*A[l,i])
    for i in range(n): # Update transformation matrix
        temp = p[i,k]
        p[i,k] = temp - s*(p[i,l] + tau*p[i,k])
        p[i,l] = p[i,l] + s*(temp - tau*p[i,l])


 maxRot = 5*(n**2) # Set limit on number of rotations
 p = identity(n)*1.0 # Initialize transformation matrix
 for i in range(maxRot): # Jacobi rotation loop
    Amax,k,l = maxElem(A)
    if Amax < tol:        
        eigvec = p.tolist()
        return diagonal(A).tolist(),eigvec
    rotate(A,p,k,l)


n = int(input("Enter the number of rows:"))
m = int(input("Enter the number of columns:"))
print("Enter the entries rowise (separated by space): ")
matrix=[]
for i in range(n):
    a=[]
    for j in range(m):
        a.append(float(input()))
    matrix.append(a)

A = np.array(matrix)
w, v = np.linalg.eig(A)
eigvec,eigval=Jacobi(A,tol = 1.0)
print(eigval)
print(eigvec)