import numpy as np
from numpy import array



seq1 = "CTAACGTTCGTC"            #   C AAT_TG A
seq2 = "AGCGTAGAATTC"           #   G AATCTG C

A = np.zeros((len(seq1)+1, len(seq2)+1),dtype=np.int16)
# scoring= array([[5, -1, 0, -1],
#                 [-1, 5, -1, 0],
#                 [0, -1, 5, -1],
#                 [-1, 0, -1, 5]])
gap = -7
scoring = [[10, -5,-5,-5],
                [-5, 10, -5, -5],
                [-5, -5, 10, -5],
                [-5, -5, -5, 10]]
nucleotides = ["A", "C", "G", "T"]


# scoring_matrix = pd.DataFrame(scoring, index = nucleotides, columns = nucleotides)



def getScore(char1,char2): #char1 ~ left, char2 ~ right 
    
    for i in range(0,len(nucleotides)):
        if char1 == nucleotides[i]:
            row = i
        if char2 == nucleotides[i]:
            col = i
    
    return scoring[row][col]        



def makeAlign(seq1,seq2):
    best = 0
    

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            # the local alignment recurrance rule:
            A[i][j] = max(A[i][j-1]   + gap,
                          A[i-1][j]   + gap,
                          A[i-1][j-1] + getScore(seq1[i-1],seq2[j-1]),
                          0)
            
            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                best_pos = []
                best_pos.append(i)
                best_pos.append(j)
    print(A)


    i = best_pos[0]
    j = best_pos[1]
    aligned_seq1 = ""
    aligned_seq2 = ""
    matchScore = 0
    
    while A[i,j] != 0:
        if A[i-1, j-1] + getScore(seq1[i-1], seq2[j-1]) == A[i,j]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            
            if seq1[i-1] == seq2[j-1]:
                matchScore += 1
            
            i -= 1
            j -= 1
            
        elif A[i-1, j] + gap == A[i,j]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "_" + aligned_seq2
            
            
            i -= 1
            
        else:
            aligned_seq1 = "_" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
#     while i > 0:
#             aligned_seq1 = "_" + aligned_seq1
#             i -= 1
#       
#     while j > 0:
#             aligned_seq2 = "_" + aligned_seq2
#             j -= 1
     
    print(i)
    print(j) 
    print(aligned_seq1)
    print(aligned_seq2)
    print(matchScore)
    return


best_pos = makeAlign(seq1, seq2)
