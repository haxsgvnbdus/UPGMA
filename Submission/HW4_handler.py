######################################################################################
#  File            : HW4_handler.py
#  Purpose         : 1> Create local alignment algorithm to maximize the score using 
#                    a scoring matrix provided
#                    2> Create a distance matrix for a set of sequences using 
#                    GLOBAL alignment.
#  Developer       : Han Ngo, Rasmi Lamichhane, Sushma Colanukudhuru 
#                    CSCI 4314 Spring 2017, HW4
##############################################################################
#   
#   Sample command line arguments to run the program: 
#   python HW4_handler.py –f dataset1.fasta –s scoring1.txt -g -1 –t dataset1.tree 
#   -f: fasta .txt input file 
#   -s: scoring matrix file
#   -g: gap penalty
#   -t: tree file which store the Newick format
#
##############################################################################
#
#   Runtime and Space Complexity:
#
#       1> Local alignment for n sequences take (nC2) comparisons, 
#       each taking up O(mxk) running time (given m and k being 
#       two sequence length and n being the number of sequences)
#
#       2> Distance matrix utilizes matrix from numpy package, which 
#       occupies O(nxn) and is used for table joining and shrinking 
#
#       3> Matrix O(nxn) is transformed to list (which is actually 
#       the lower triangular matrix) to optimize running time and 
#       program memory. Instead of nxn, the list stores neccessary 
#       nxn/2 elements and run O(nxn) at most for adjusting the table as
#
##############################################################################

import sys, getopt, math, os.path
from distanceMatrix import *
from formNewick import *

##############################################################################
# 
#   Specification of the command line to run the program
#               
##############################################################################


def usage():
    print "An input FASTA file (-f)"

##############################################################################
# 
#   Main function in this handler .py file takes in argument and manipulate 
#   them accordantly by importing all functions in related files
#               
##############################################################################
  
def main(argv):    
    gap = 0
    
    try:
        opts, args = getopt.getopt(argv, "hf:s:g:t:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        if opt == '-f':
            if os.path.isfile(arg):
                print "\nFile: " + arg
                sequences = readInput(arg)    
        if opt =='-g':
            gap = int(arg)
        if opt == '-s':
            if os.path.isfile(arg):
                matrixFile=arg
                
                #Read scoring matrix file and storing content 
                scores, nus = readmatrix(matrixFile)    
                    
                #Create distance matrix     
                distMatrix = makeDistanceMatrix(sequences, nus, scores, gap)
                
                #convert Matrix to list (lower triangular matrix) 
                listDistance = matrixToList(distMatrix)
                
                #shorten sequence names
                seq = naming(sequences)
                
                #generate Newick format 
                result = createNewick(listDistance, seq)
            
            # If the input file given does not exist, print an error message and exit the program
            else:
                print "ERROR! Input file must exist. Correct usage:"
                usage()
                sys.exit(2)
                
    if opt == '-t':
        if os.path.isfile(arg):
            treeFile = arg
            f = open(treeFile, 'w')
            f.write(result)
        
        # If the input file given does not exist, exit the program
        else:
                print "ERROR! Tree file must exist. Correct usage:"
                sys.exit(2)

##############################################################################
# 
#   Store all sequences in a list
#   *Sequence names are neglected. Note the case that long sequences are
#   two lines. 
#               
##############################################################################

def readInput(inFile):
    f = open(inFile, 'r')
    lines =  f.read().splitlines()
    temp_seq = ""
    sequences = []    

    for i in range(len(lines)):    
        if lines[i].startswith(">"):
        temp_seq = ""        #renew the temporary sequence to store
        else:
        temp_seq = temp_seq + lines[i]
        try:    
            # check that if new sequence starts the next line, append the last one    
            if lines[i+1].startswith(">"):
                sequences.append(temp_seq)
        except:
            sequences.append(temp_seq)
            
    return sequences


##############################################################################
# 
#   Reading scoring matrix file and return score values as a matrix and 
#   its accodant nucleotides "G", "C", "T", "A" as the first element
#   of any row.
#               
##############################################################################


def readmatrix(inFile):
    
    nus = []
    with open(inFile) as textFile:
        scores = [line.split() for line in textFile]
    for row in scores:
        nus.append(row[0])
        del row[0]

    return scores, nus

##############################################################################
# 
#   Trigger the program
#               
##############################################################################

if __name__ == "__main__":
    main(sys.argv[1:])