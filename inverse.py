import sys

def inverse(A):
    def print_matrix(Title, M):
        print(Title)
        for row in M:
            print([round(x, 3)+0 for x in row])


    def print_matrices(Action, Title1, M1, Title2, M2):
        print(Action)
        print(Title1, '\t'*int(len(M1)/2)+"\t"*len(M1), Title2)
        for i in range(len(M1)):
            row1 = ['{0:+7.3f}'.format(x) for x in M1[i]]
            row2 = ['{0:+7.3f}'.format(x) for x in M2[i]]
            print(row1, '\t', row2)


    def zeros_matrix(rows, cols):
        A = []
        for i in range(rows):
            A.append([])
            for j in range(cols):
                A[-1].append(0.0)

        return A


    def copy_matrix(M):
        rows = len(M)
        cols = len(M[0])

        MC = zeros_matrix(rows, cols)

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

        C = zeros_matrix(rowsA, colsB)

        for i in range(rowsA):
            for j in range(colsB):
                total = 0
                for ii in range(colsA):
                    total += A[i][ii] * B[ii][j]
                C[i][j] = total
        return C


    # A = [[6, 1, -1],[0, 7, 0],[3, -1, 2]]
    # I=[[1,0,0],[0,1,0],[0,0,1]]



    def identity_matrix(x):
        matrix = [[0 for j in range(x)] for i in range(x)]
        for i in range(x):
            matrix[i][i] = 1
        return matrix

    I=identity_matrix(len(A))
    # print_matrices('','A Matrix', A, 'I Matrix', I)

    AM = copy_matrix(A)
    IM = copy_matrix(I)
    n = len(AM)

    exString = """
    Since the matrices won't be the original A and I as we start row operations, 
        the matrices will be called: AM for "A Morphing", and IM for "I Morphing" 
    """


    fd = 0 
    fdScaler = 1. / AM[fd][fd]
    print(AM)
    for j in range(n):
        AM[fd][j] = fdScaler * AM[fd][j]
        IM[fd][j] = fdScaler * IM[fd][j]
        print(AM)
        print(IM)
        


    n = len(A)
    indices = list(range(n))

    for i in indices[0:fd] + indices[fd+1:]:
        crScaler = AM[i][fd]
        for j in range(n):
            AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
            IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
        

    indices = list(range(n))
    for fd in range(1,n):
        fdScaler = 1.0 / AM[fd][fd]
        for j in range(n):
            AM[fd][j] *= fdScaler
            IM[fd][j] *= fdScaler
        
        string1 = '\nUsing the matrices above, Scale row-{} of AM and IM by '
        string2 = 'diagonal element {} of AM, which is 1/{:+.3f}.\n'
        stringsum = string1 + string2
        val1 = fd+1; val2 = fd+1
        Action = stringsum.format(val1,val2,round(1./fdScaler,3))
        
        
        for i in indices[:fd] + indices[fd+1:]:
            crScaler = AM[i][fd]
            for j in range(n):
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
                print(AM)
            
            string1 = 'Using the matrices above, subtract {:+.3f} * row-{} of AM from row-{} of AM, and \n'
            string2 = '\tsubtract {:+.3f} * row-{} of IM from row-{} of IM\n'
            val1 = i+1; val2 = fd+1
            stringsum = string1 + string2
            Action = stringsum.format(crScaler, val2, val1, crScaler, val2, val1)
            


    print_matrix('Inverse of input matrix: ', matrix_multiply(IM,I))
    x=matrix_multiply(IM,I)
    return x

    # print_matrix('Proof of Inversion', matrix_multiply(A,IM))

    
n = int(input("Enter the number of rows:"))
m = int(input("Enter the number of columns:"))
print("Enter the entries rowise (separated by space): ")
A=[]
for i in range(n):
    a=[]
    for j in range(m):
        a.append(int(input()))
    A.append(a)

inv=inverse(A)
print(inv)
