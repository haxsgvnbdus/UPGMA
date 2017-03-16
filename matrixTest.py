'''
Created on Mar 15, 2017

@author: Hannie
'''

def matrixToList(table):
    matrix = []
    for i in range(0, len(A[0])):
        row = []
        for j in range(0, len(A[0])):
            if i>j:
                row.append(A[i][j])
        matrix.append(row)
#     print(pos)
    print(matrix)
    
    
    