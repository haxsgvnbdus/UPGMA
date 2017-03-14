'''
Created on Mar 11, 2017

@author: Hannie
'''

import re

#Part 2: UPGMA | Supposed given distance matrix 
#Sample in textbook

def findMinCell(table):
    # Set default to infinity
    min_cell = float("inf")
    x, y = -1, -1

    # Go through every cell, looking for the lowest
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j

    # Return the x, y co-ordinate of cell
#     print("Minimum value: " +  str(min_cell) + " at " + str(x) + " " + str(y))
    return x, y, min_cell




def joinSequences(sequences, a, b, distance):      #   Merge 2 sequences 
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a
    distance = round(distance,1)
    
    try: 
        if (sequences[a][len(sequences[a])-1:]) == ")":
            print(sequences[a][len(sequences[a])-1:])
            result = re.findall("\d+\.?\d*", sequences[a])
            distance_to_base = float(result[len(result)-1])
            print(distance_to_base)
            finalDist = round(distance/2 - distance_to_base,2)
            print("distance: " + str(distance))
            print("final distance: " + str(finalDist))
            
        
        # Join the sequences in the first index
            sequences[a] = "(" + sequences[a] + ":" + str(finalDist) + "," + sequences[b] + ":" +  str(distance/2) + ")"
            print(sequences[a])
            del sequences[b]
            return

    except:  
        pass
    sequences[a] = "(" + sequences[a] + ":" + str(distance/2) + "," + sequences[b] + ":" +  str(distance/2) + ")"
    print(sequences[a])
    
    del sequences[b]

    # Remove the (now redundant) label in the second index
    

#   Joins the entries of a table on the cell (a, b) by averaging their data entries
def shrinkTable(table, a, b):
    # Swap if the indices are not ordered
    if b < a:
#         print("swap " + str(a) + " " + str(b))
        a, b = b, a
    
    # For the lower index, reconstruct the entire row (a, i), where i < A 
    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i])/2)
    table[a] = row

    #Lower matrix: row index > column index    
    # Reconstruct the entire column (i, A), where i > A
    #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
    for i in range(a+1, b):
        table[i][a] = (table[i][a]+table[b][i])/2
        
    #   Finish up the rest of the values from row i
    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2
        
        
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]




def UPGMA(table, sequences):
    
    while len(sequences) > 1:
        #min_cell is also twice times the distance from the base to its common ancestor
        x, y, min_cell = findMinCell(table)
#         print("x: " + str(x) + " y: " + str(y) + "distance: " + str(min_cell/2))
        
        shrinkTable(table, x, y)
        joinSequences(sequences, x, y, min_cell)

    return sequences[0]


def testingSample(start, end):
    sequences = []
    for i in range(ord(start), ord(end)+1):
        sequences.append(chr(i))
    print(sequences)
    return sequences


M_sequences = testingSample("A", "G")   
M = [
    [],                         #A
    [19],                       #B
    [27, 31],                   #C
    [8, 18, 26],                #D
    [33, 36, 41, 31],           #E
    [18, 1, 32, 17, 35],        #F
    [13, 13, 29, 14, 28, 12]    #G
    ]


UPGMA(M, M_sequences)  #should output: '((((A,D),((B,F),G)),C),E)'