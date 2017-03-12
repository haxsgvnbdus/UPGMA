'''
Created on Mar 11, 2017

@author: Hannie
'''

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
    print("Minimum value: " +  str(min_cell) + " at " + str(x) + " " + str(y))
    return x, y



def joinSequences(sequences, a, b):      #   Merge 2 sequences 
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the sequences in the first index
    sequences[a] = "(" + sequences[a] + "," + sequences[b] + ")"
    print(sequences[a])
    # Remove the (now redundant) label in the second index
    del sequences[b]


#Joins the entries of a table on the cell (a, b) by averaging their data entries
def shrinkTable(table, a, b):



def UPGMA(table, sequences):

    while len(sequences) > 1:
    
        x, y = findMinCell(table)
        shrinkTable(table, x, y)
        joinSequences(sequences, x, y)

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