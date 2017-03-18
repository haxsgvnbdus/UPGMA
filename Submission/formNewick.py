################################################################################
#  File            : formNewick.py
#  Purpose         : 1> Use lower triangle matrix to do UPGMA neighbor joining
#                    2> Create trees via www.trex.uqam.ca/index.php?action=newick
################################################################################

import re

##############################################################################
# 
#   Find the smallest element in concurrent table
#               
##############################################################################

def findMinCell(table):
    # Set default to min value and its co-ordinates (x,y)
    min_cell = float("inf")
    x, y = -1, -1

    # Go through every cell, looking for the lowest
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j

    return x, y, min_cell


##############################################################################
# 
#   Merging 2 sequences and calculating updated distance to the common 
#   ancestor. Specifically the local variable will store the subformat of 
#   Newick and keep updated till we get to the end
#               
##############################################################################

def joinSequences(seqs, x, y, dist):      
    # Swap if the indices are not ordered
    if y < x:
        x, y = y, x

    updateX = updateY = 0
    try: 
        updateX = updateDist(seqs[x], dist) 
        updateY = updateDist(seqs[y], dist) 
        
    except:  
        pass
    
    # Join the sequences in the first index
    seqs[x] = "(" + seqs[x] + ":" + str(updateX) + "," +     \
                    seqs[y] + ":" +  str(updateY) + ")"
    
    # Remove the (now redundant) label in the second index
    del seqs[y]          
    return seqs[x] 

##############################################################################
# 
#   Joins the entries of a table on the cell (x, y) by averaging 
#   their data entries and delete the redundant (already joined) ones
#               
##############################################################################
def shrinkTable(distance_list, x, y):
    # Swap if the indices are not ordered
    if y < x:
        x, y = y, x
    
    # For the lower index, reconstruct the entire row (x, i), where i < x 
    row = []
    for i in range(0, x):
        row.append((distance_list[x][i] + distance_list[y][i])/2)
    distance_list[x] = row

    #Lower matrix: row index > column index    
    
    # Reconstruct the entire column (i, A), where i > A
    #Note: Since the matrix is lower triangular, row y only contains values for indices < y
    for i in range(x+1, y):
        distance_list[i][x] = (distance_list[i][x]+distance_list[y][i])/2
        
    #   Finish up the rest of the values from row i
    for i in range(y+1, len(distance_list)):
        distance_list[i][x] = (distance_list[i][x]+distance_list[i][y])/2
        
        # Remove the (now redundant) second index column entry
        del distance_list[i][y]

    # Remove the (now redundant) second index row
    del distance_list[y]

##############################################################################
# 
#   The main caller in this helper of creating Newick format
#               
##############################################################################
def createNewick(distance_list, seqName):

    while len(seqName) > 1:
        
        #Note: min_cell is 2x distance from base to its common ancestor
        x, y, min_cell = findMinCell(distance_list)
        
        shrinkTable(distance_list, x, y)
        result = joinSequences(seqName, x, y, min_cell)

    print(result)
    return result

##############################################################################
# 
#   Updating distance to common ancestor 
#               
##############################################################################
def updateDist(seq,dist):
    
    list_dist = re.findall("\d+\.?\d*", seq)

    #There are subgroups
    if (seq[len(seq)-1:]) == ")":
        list_dist = sorted(list_dist, key = float, reverse = True)
        new_dist = round(dist/2 - float(list_dist[0]),5)
        return new_dist
    
    #There is no subgroup
    else: 
        return round(dist/2,5)

##############################################################################
# 
#    As said in HW4_handler.py, this is to convert matrix to list for 
#    easy reference, convenient cell deletion, and optimize space
#               
##############################################################################
def matrixToList(A):
    list = []
    for i in range(0, len(A[0])):
        row = []
        for j in range(0, len(A[0])):
            if i>j:
                row.append(A[i][j])
        list.append(row)
    
    return list

##############################################################################
# 
#    Keep the first 8 characters as ID for the sequences to put in dendrogram
#
##############################################################################    
def naming(sequences):
    renamed = []
    for i in range(0, len(sequences)):
        renamed.append(sequences[i][0:8])
    
    return renamed
