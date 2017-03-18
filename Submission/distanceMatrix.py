################################################################################
#  File            : distanceMatrix.py
#  Purpose         : 1> Implement local alignment for every two sequences
#                    2> Create distance matrix 
################################################################################

import numpy as np
from numpy import array

##############################################################################
# 
#    Create distance matrix (1 vs 2, 1 vs 3, etc, 2 vs 3, etc, N vs N) 
#    using a local alignment with customized scoring matrix.
#    Print distance matrix (all values are between 0-1, use 5 decimal places)
#    meaning (1 - similarity = 1 - #matches/length of sequence)
#               
##############################################################################
def makeDistanceMatrix(sequences, nus, scores, gap): 
    distMatrix = np.zeros((len(sequences), len(sequences)))

    #Set 5 decimals for the division result
    np.set_printoptions(formatter={'float': '{: 0.5f}'.format})

    for i in range(len(sequences)):
        for j in range(len(sequences)):
            #comparing itself, the distance scores 0 
            if i == j:    
                distMatrix[i][j] = 0
                
            #comparing with every other sequence, score 1 - similarity as said in header
            elif i<j:     
                distMatrix[i][j] = 1- makeAlign(sequences[i] ,sequences[j], nus, scores, gap)
                distMatrix[j][i] = distMatrix[i][j]
                
    print(distMatrix)
    return distMatrix


##############################################################################
#     Applying Dynamic Programming and local alignment algorithm to 
#     fill in all table cells. 
#     matchScore is to calculate every nucleotide that match
#     gapScore is to count how many gap are there in the alignment 
#     There two are used to calculate the distance matrix 
#
##############################################################################

def makeAlign(seq1,seq2, nus, scores, gap):
    best = 0
    
    #Initiation
    A = np.zeros((len(seq1)+1, len(seq2)+1),dtype=np.int16)


    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            
            # the local alignment recurrance rule:
            A[i][j] = max(A[i][j-1]   + gap,
                          A[i-1][j]   + gap,
                          A[i-1][j-1] + int(getScore(nus, scores, seq1[i-1],seq2[j-1])),
                          0)
            
            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                best_pos = []
                best_pos.append(i)
                best_pos.append(j)
    #print(A)

    #backtracking part
    i = best_pos[0]
    j = best_pos[1]
    matchScore = 0
    gapScore = 0
    
    while A[i,j] != 0:
        
        #Priority: check from diagonal
        if A[i-1, j-1] + int(getScore(nus, scores, seq1[i-1], seq2[j-1])) == A[i,j]:            
            if seq1[i-1] == seq2[j-1]:
                matchScore += 1
            
            i -= 1
            j -= 1
            
        #Check from top
        elif A[i-1, j] + gap == A[i,j]:
            gapScore += 1
            i -= 1
        
        #Check from left    
        else:
            gapScore += 1
        j -= 1
    

    return matchScore/float(len(seq1)+gapScore)
 
##############################################################################
# 
#    getScore return values as whether the two nucleotides are the same or not
#               
##############################################################################
def getScore(nus, scores, char1, char2): 
    #Initiate row and column index as zero
    r = 0
    c = 0

    for i in range(0,len(nus)):
        if char1 == nus[i]:
            r = i
        if char2 == nus[i]:
            c = i

    return scores[r][c]

