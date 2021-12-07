# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 14:15:27 2021

@author: Matthew Chan
"""

##############################################################################
# OBTAIN SEQEUNCES
##############################################################################

import tkinter as tk
import tkinter.filedialog as fd

root = tk.Tk()
files = fd.askopenfilenames(parent=root, title='Choose a file')
root.mainloop()

DNA = open(files[0])
DNA_name = DNA.readline()
DNA_seq = DNA.read()
DNA.close()

prot = open(files[1])
prot_name = prot.readline()
prot_seq = prot.read()
prot.close()

DNA_seq = DNA_seq.replace("\n", "")
prot_seq = prot_seq.replace("\n", "")

##############################################################################
# GET TRANSLATED DNA to PROTEIN SEQUENCE
##############################################################################

def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    import re
    
    indexes = [m.start() for m in re.finditer('ATG', seq)] #Finds all ATG
    prot_list = []
    #Loop through DNA sequence to match amino acid to codons
    for index in indexes:
        protein = ""
        for i in range(index, len(seq), 3):
            codon = seq[i:i+3]
            protein += table[codon]
            if codon == "TAA" or codon == "TAG" or codon == "TGA":
                break
        if protein[-1] == "_":
            prot_list.append(protein[0:len(protein) - 1])
        
    long_prot = max(prot_list, key=len)
    
    return long_prot

DNA_prot = translate(DNA_seq) #Translate DNA sequence to protein

##############################################################################
# DUMB STUFF TO GET BLOSUM MATRIX
##############################################################################

Blosum = open("BLOSUM62.txt")
for x in range(0, 6):
    Blosum.readline()

BLOSUM62_matrix = []
BLOSUM62_names = []

#Obtains the animo acide order regarding their letter name
temp = Blosum.readline()
temp = temp.replace(" ", "")
for y in temp:
    BLOSUM62_names.append(y)
BLOSUM62_names = BLOSUM62_names[0:len(BLOSUM62_names) - 1]

#Creates a matrix containing values for BLOSUM62
for z in range(0, len(BLOSUM62_names)):
    line_array = []
    line = Blosum.readline()
    line = line.replace(" ", "")
    line = line[1:len(line) - 1]
    count = 0
    while count < len(line):
        if line[count] == "-":
            line_array.append(int(line[count:count+2]))
            count += 2
        else:
            line_array.append(int(line[count]))
            count += 1
    BLOSUM62_matrix.append(line_array)
    
Blosum.close()

##############################################################################
# CREATING SCORE MATRIX
##############################################################################

gap_o = -10 #gap opening penalty
gap_e = -2 #gap extension penalty
score_matrix = []
max_index = [0, 0, 0] # x, y, value

#Fills in first row of score matrix with zeroes because local alignment
first_row = []
for x in range(0, len(prot_seq) + 1):
    first_row.append(0)
score_matrix.append(first_row)

#Initializes matrix to be same size and adds 0 for first column while the rest
#are blank
for y in range(0, len(DNA_prot)):
    temp = [0]
    for x in range(0, len(prot_seq)):
        temp.append("")
    score_matrix.append(temp)

#Initializes the prev variable which is used to determine whether previous
#score came from gap or match
if BLOSUM62_names.index(prot_seq[0]) == BLOSUM62_names.index(DNA_prot[0]):
    prev = 0
else:
    prev = 1

#Fills score matrix depending on highest value from gap or BLOSUM62 match
for y in range(0, len(prot_seq)):
    for x in range(0, len(DNA_prot)):
        aa_x = prot_seq[y] #Amino acid from protein sequence
        aa_y = DNA_prot[x] #Amino acid from translated DNA sequence
        
        #Finds the BLOSUM62 score of the two proteins at certain position
        index_x = BLOSUM62_names.index(aa_x)
        index_y = BLOSUM62_names.index(aa_y)
        
        match = BLOSUM62_matrix[index_x][index_y]
        
        #Determines all possible scores depending on whether it's a gap
        #opening or extension
        if prev == 0:
            val1 = score_matrix[x][y + 1] + gap_o
            val2 = score_matrix[x + 1][y] + gap_o
        else:
            val1 = score_matrix[x][y + 1] + gap_e
            val2 = score_matrix[x + 1][y] + gap_e
        val3 = score_matrix[x][y] + match
        
        #Determines max value of all possible scores
        max_val = max(val1, val2, val3)
        
        #Determines prev variable whether match came from match or gap
        if max_val == val3:
            prev = 0
        else:
            prev = 1
        
        #Fills score matrix with the max value
        if max_val < 0:
            score_matrix[x+1][y+1] = 0
        else:
            score_matrix[x+1][y+1] = max_val
        
        #Keeps track of position and max value for traceback purposes
        if max_val >= max_index[2]:
            max_index[0] = x + 1
            max_index[1] = y + 1
            max_index[2] = max_val 

##############################################################################
# TRACEBACK
##############################################################################

pos = [max_index[0],max_index[1]] #position variable

#Start of the alignment with max value amino acid match
DNA_align = DNA_prot[pos[0]-1]
prot_align = prot_seq[pos[1]-1]

max_num = max_index[2] #Max value

#Traceback when max value is greater than zero
while max_num > 0:
    
    #Determines the three possible values for traceback
    num1 = score_matrix[pos[0]-1][pos[1]-1]
    num2 = score_matrix[pos[0]][pos[1]-1]   #up?
    num3 = score_matrix[pos[0]-1][pos[1]]   #left?
    
    #Adds an astrix to denote the path
    score_matrix[pos[0]][pos[1]] = "*" + str(score_matrix[pos[0]][pos[1]])
    
    #If all values are 0 it has reached end of the traceback
    if num1 == 0 and num2 == 0 and num3 == 0:
        break
    
    #When the highest value is diagonal it prioritizes it
    if num1 >= num2 and num1 >= num3:
        pos = [pos[0]-1, pos[1]-1]
        
        DNA_align = DNA_prot[pos[0]-1] + DNA_align
        prot_align = prot_seq[pos[1]-1] + prot_align
        
        max_num = num1
    
    #Highest match is the upper value
    elif num2 >= num3:
        pos = [pos[0], pos[1]-1]
        
        DNA_align = "-" + DNA_align
        prot_align = prot_seq[pos[1]-1] + prot_align
        
        max_num = num2
    
    #Highest match is the left value
    else:
        pos = [pos[0]-1, pos[1]]
        
        DNA_align = DNA_prot[pos[0]-1] + DNA_align
        prot_align = "-" + prot_align
        
        max_num = num3

##############################################################################
# OUTPUT SCORE MATRIX AND SEQUENCE ALIGNMENT
##############################################################################

#Creation of csv file of the score matrix
file = open("Local_ScoreMatrix.csv", "w")

file.write(", ,")
for x in range(0, len(prot_seq)):
    file.write(prot_seq[x] + ",")
file.write("\n")

DNA_prot = " " + DNA_prot

for i in range(0, len(score_matrix)):
    file.write((DNA_prot[i]) + ",")
    for j in range(0, len(score_matrix[0])):
        file.write(str(score_matrix[i][j]) + ",")
    file.write("\n")
file.close()

#Adds the position in the sequence to the alignment 
DNA_align_temp = DNA_align.replace("-", "")
prot_align_temp = prot_align.replace("-", "")
DNA_start = DNA_prot.index(DNA_align_temp)
prot_start = prot_seq.index(prot_align_temp)
if len(str(prot_start)) > len(str(DNA_start)):
    spaces = " " * (len(str(prot_start)) - len(str(DNA_start)))
    spaces2 = " " * (len(str(prot_start + len(prot_align_temp))) - len(str(DNA_start + len(DNA_align_temp))))
    DNA_align = str(DNA_start) + spaces + "   " + DNA_align + "   " + spaces2 + str(DNA_start + len(DNA_align_temp) - 1)
    prot_align = str(prot_start + 1) + "   " + prot_align + "   " + str(prot_start + len(prot_align_temp))
else:
    spaces = " " * (len(str(DNA_start)) - len(str(prot_start)))
    spaces2 = " " * (len(str(DNA_start + len(DNA_align_temp))) - len(str(prot_start + len(prot_align_temp))))
    prot_align = str(prot_start + 1) + spaces + "   " + prot_align + "   " + spaces2 + str(prot_start + len(prot_align_temp))
    DNA_align = str(DNA_start) + "   " + DNA_align + "   " + str(DNA_start + len(DNA_align_temp) - 1)

#Creates file with some info and the alignment
file2 = open("Local_SequenceInfo.txt", "w")

file2.write("DNA:\t\t" + DNA_name[1:len(DNA_name)])
file2.write("Protein:\t" + prot_name[1:len(prot_name)])
file2.write("\nTranslated DNA: " + DNA_prot)
file2.write("\nProtein: " + prot_seq)
file2.write("\nUsed BLOSUM62 for substitution matrix\n")
file2.write("Gap Opening = " + str(gap_o))
file2.write("\nGap Extension = " + str(gap_e))
file2.write("\n\n\nALIGNMENT\n\n")
file2.write("DNA:\t\t" + DNA_align + "\n\t\t")
for k in range(0, len(DNA_align)):
    if DNA_align[k] == "-" or prot_align[k] == "-" or DNA_align[k] != prot_align[k] or DNA_align[k] == " " or DNA_align[k].isnumeric() is True:
        file2.write(" ")
    else:
        file2.write("|")
file2.write("\nProtein:\t" + prot_align + "\n")
file2.write("\nNote: Look at Local_ScoreMatrix.csv for full scoring matrix")

file2.close()





    
    
    
    
    
    