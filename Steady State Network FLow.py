from math import sqrt
import numpy as np
import sys

##########################################################################################################
def Jacobi(A,tol = 1.0e-9): 
    # Jacobi method
    # 'A' Matrix is numpy Matrix formatted in Jacobi function
    #  To find the largest element of off-diagonal A[d,e]
    def Max_ELement(A):
        n = len(A)
        A_Maximum = 0.0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(A[i,j]) >= A_Maximum:
                    A_Maximum = abs(A[i,j])
                    k = i; l = j
        return A_Maximum,k,l

    # Function Rotate_Matrix Matrix makes A[d,e]=0 and defines the Rotation Matrix
    def Rotate_Matrix(A,r,d,e):
        n = len(A)
        A_Difference = A[e,e] - A[d,d]
        if abs(A[d,e]) < abs(A_Difference)*1.0e-36: t = A[d,e]/A_Difference
        else:
            phi = A_Difference/(2.0*A[d,e])
            t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0/sqrt(t**2 + 1.0); s = t*c
        Tau = s/(1.0 + c)
        Temp = A[d,e]
        A[d,e] = 0.0
        A[d,d] = A[d,d] - t*Temp
        A[e,e] = A[e,e] + t*Temp
        for i in range(d): # This is the case of i < d
            Temp = A[i,d]
            A[i,d] = Temp - s*(A[i,e] + Tau*Temp)
            A[i,e] = A[i,e] + s*(Temp - Tau*A[i,e])
        for i in range(d+1,e): # This is the case of d < i < e
            Temp = A[d,i]
            A[d,i] = Temp - s*(A[i,e] + Tau*A[d,i])
            A[i,e] = A[i,e] + s*(Temp - Tau*A[i,e])
        for i in range(e+1,n): # This is the case of i > e
            Temp = A[d,i]
            A[d,i] = Temp - s*(A[e,i] + Tau*Temp)
            A[e,i] = A[e,i] + s*(Temp - Tau*A[e,i])
        for i in range(n): # This 'for' loop updates transformation matrix
            Temp = r[i,d]
            r[i,d] = Temp - s*(r[i,e] + Tau*r[i,d])
            r[i,e] = r[i,e] + s*(Temp - Tau*r[i,e])

    Max_Rot = 5*(n**2) # Setting the limit on number of rotations
    p = np.identity(n)*1.0 # Initializing the transformation matrix
    for i in range(Max_Rot): # Jacobi Rotation 'for' loop
        A_Maximum,k,l = Max_ELement(A)
        if A_Maximum < tol:        
            eigvec = p.tolist()
            return eigvec,np.diagonal(A).tolist()
        Rotate_Matrix(A,p,k,l)

################################################################################################################################################

def Inverse(A): 
    # Inverse function contains a list 'A'
    def Zero_Matrix(rows, cols):
        A = []
        for i in range(rows):
            A.append([])
            for j in range(cols):
                A[-1].append(0.0)

        return A


    def Copy_Matrix(M):
        rows = len(M)
        cols = len(M[0])

        MC = Zero_Matrix(rows, cols)

        for i in range(rows):
            for j in range(rows):
                MC[i][j] = M[i][j]

        return MC

    def matrix_multiply(A, B):
        rowsA = len(A)
        colsA = len(A[0])

        rowsB = len(B)
        colsB = len(B[0])

        if colsA != rowsB:
            print('Number of A columns must equal number of B rows.')
            sys.exit()

        C = Zero_Matrix(rowsA, colsB)

        for i in range(rowsA):
            for j in range(colsB):
                total = 0
                for ii in range(colsA):
                    total += A[i][ii] * B[ii][j]
                C[i][j] = total

        return C


    def Identity_Matrix(x):
        matrix = [[0 for j in range(x)] for i in range(x)]
        for i in range(x):
            matrix[i][i] = 1
        return matrix

    I=Identity_Matrix(len(A))
    AM = Copy_Matrix(A)
    IM = Copy_Matrix(I)
    n = len(AM)

    # A and I matrix are not original matrix as we start the row operations.
    # So, the matrices will be called AM for 'A Morphing' and 'IM for I Morphing'
    fd = 0 
    fdScaler = 1. / AM[fd][fd]

    for j in range(n):
        AM[fd][j] = fdScaler * AM[fd][j]
        IM[fd][j] = fdScaler * IM[fd][j]

    n = len(A)
    Indices = list(range(n))

    for i in Indices[0:fd] + Indices[fd+1:]:
        crScaler = AM[i][fd]
        for j in range(n):
            AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
            IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
        

    Indices = list(range(n))
    for fd in range(1,n):
        fdScaler = 1.0 / AM[fd][fd]
        for j in range(n):
            AM[fd][j] *= fdScaler
            IM[fd][j] *= fdScaler
        
        
        for i in Indices[:fd] + Indices[fd+1:]:
            crScaler = AM[i][fd]
            for j in range(n):
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
    x=matrix_multiply(IM,I)
    return x

###############################################################################################################################

def diagonilization(eigen_val,eigen_vec,iter):
    eigvec_Inverse=Inverse(eigen_vec)
    def matrix2(y):
        arr = [[0 for j in range(1)] for i in range(y)]
        x=0
        for i in range(y):
            arr[i][0]=1/y
            x+=1
            if x>y:
                break
        return arr

    def diagonal_eigval(eigen_val):
        arr = [[0 for j in range(len(eigen_val))] for i in range(len(eigen_val))]
        x=0
        for i in eigen_val:
            arr[x][x]=pow(i,iter)
            x+=1
            if x>len(eigen_val):
                break
        return arr

    def Zero_Matrix(rows, cols):
        A = []
        for i in range(rows):
            A.append([])
            for j in range(cols):
                A[-1].append(0.0)

        return A

    def matrix_multiply(A, B):
        rowsA = len(A)
        colsA = len(A[0])

        rowsB = len(B)
        colsB = len(B[0])

        if colsA != rowsB:
            print('Number of A columns must equal number of B rows.')
            sys.exit()

        C = Zero_Matrix(rowsA, colsB)

        for i in range(rowsA):
            for j in range(colsB):
                total = 0
                for ii in range(colsA):
                    total += A[i][ii] * B[ii][j]
                C[i][j] = total
        return C


    dia_eigval=diagonal_eigval(eigen_val)
    x=matrix_multiply(eigen_vec,dia_eigval)
    z=matrix_multiply(eigvec_Inverse,matrix2(len(eigvec_Inverse)))
    y=matrix_multiply(x,z)

    return y
#################################################################################################################################################
# Input Begins Here: 
print("Instruction: ")
print("1) The Input matrix of the network should be 'Symmetric Matrix'.")
print("2) The sum of Column elments of the Matrix should be 'Unity'.")
print("3) Input Non-Negative values.")
n = int(input("Enter the number of rows:"))
m = int(input("Enter the number of columns:"))
print("Enter the entries rowise(one element at a time): ")
matrix=[]
for i in range(n):
    a=[]
    for j in range(m):
        a.append(float(input()))
    matrix.append(a)

A_Matrix=np.array(matrix)
A_list=matrix
print("Input Matrix: ")
for i in A_list:
    print(i)
print("")
eigen_vec,eigen_val=Jacobi(A_Matrix,tol = 1.0e-4)
print("Eigen Values of the Matrix: ")
print(eigen_val)
print("")
print("Eigen Vectors of the Matrix(Column wise): ")
for i in eigen_vec:
    print(i)
print("")
eigvec_Inverse=Inverse(eigen_vec)
print("Inverse of Eigen Vector Matrix(Column wise):")
for i in eigvec_Inverse:
    print(i)
print("")
dia1=diagonilization(eigen_val,eigen_vec,100)
print("Steady State of Network Flow (after 100 iterations): ")
for i in dia1:
    print(i)
print("")
dia2=diagonilization(eigen_val,eigen_vec,200)
print("Steady State of Network Flow (after 200 iterations): ")
for i in dia2:
    print(i)