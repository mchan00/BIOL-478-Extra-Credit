# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:07:30 2021

@author: Matthew Chan
"""
'''
Look at LocalAlign_code for more comments since it's esentially the same
'''
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
    
    indexes = [m.start() for m in re.finditer('ATG', seq)]
    prot_list = []
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

DNA_prot = translate(DNA_seq)

##############################################################################
# DUMB STUFF TO GET BLOSUM MATRIX
##############################################################################

Blosum = open("BLOSUM62.txt")
for x in range(0, 6):
    Blosum.readline()

BLOSUM62_matrix = []
BLOSUM62_names = []

temp = Blosum.readline()
temp = temp.replace(" ", "")
for y in temp:
    BLOSUM62_names.append(y)
BLOSUM62_names = BLOSUM62_names[0:len(BLOSUM62_names) - 1]

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

gap_o = -10
gap_e = -2
score_matrix = []
max_index = [0, 0, 0] # x, y, value

first_row = []
for x in range(0, len(prot_seq) + 1):
    first_row.append(0)
score_matrix.append(first_row)

for y in range(0, len(DNA_prot)):
    temp = [0]
    for x in range(0, len(prot_seq)):
        temp.append("")
    score_matrix.append(temp)

if BLOSUM62_names.index(prot_seq[0]) == BLOSUM62_names.index(DNA_prot[0]):
    prev = 0
else:
    prev = 1

for y in range(0, len(prot_seq)):
    for x in range(0, len(DNA_prot)):
        aa_x = prot_seq[y]
        aa_y = DNA_prot[x]
        
        index_x = BLOSUM62_names.index(aa_x)
        index_y = BLOSUM62_names.index(aa_y)
        
        match = BLOSUM62_matrix[index_x][index_y]
        
        if prev == 0:
            val1 = score_matrix[x][y + 1] + gap_o
            val2 = score_matrix[x + 1][y] + gap_o
        else:
            val1 = score_matrix[x][y + 1] + gap_e
            val2 = score_matrix[x + 1][y] + gap_e
        val3 = score_matrix[x][y] + match
        
        max_val = max(val1, val2, val3)
        
        if max_val == val3:
            prev = 0
        else:
            prev = 1
        
        if max_val < 0:
            score_matrix[x+1][y+1] = 0
        else:
            score_matrix[x+1][y+1] = max_val

##############################################################################
# TRACEBACK
##############################################################################

import copy

num_align = 5 #Number of desired alignments

align_DNA_all = [] #Holds all the DNA alignments
align_prot_all = [] #Holds all the protein alignments

#Find maximum value in a matrix to determine position and max value in score 
#matrix
def findMax(newMatrix):
    maxval = [0,0,0]
    for x in range(0, len(DNA_prot)):
        for y in range(0, len(prot_seq)):
            value = newMatrix[x][y]
            if value >= maxval[2]:
                maxval[2] = value
                maxval[0] = x + 1
                maxval[1] = y + 1
                
    return maxval

copy_matrix = copy.deepcopy(score_matrix) #Make a copy of the score matrix

#Iterate through number of alignments desired
for i in range(0, num_align):
    max_index = findMax(score_matrix)
    
    pos = [max_index[0],max_index[1]]
    
    DNA_align = DNA_prot[pos[0]-1]
    prot_align = prot_seq[pos[1]-1]
    
    max_num = max_index[2]
    
    while max_num > 0:
        
        num1 = score_matrix[pos[0]-1][pos[1]-1]
        num2 = score_matrix[pos[0]][pos[1]-1]   #up?
        num3 = score_matrix[pos[0]-1][pos[1]]   #left?
        
        copy_matrix[pos[0]][pos[1]] = "(" + str(i+1) + ")" + str(copy_matrix[pos[0]][pos[1]])
        score_matrix[pos[0]][pos[1]] = 0 #Make value 0 so it will not be repeated
        
        if num1 == 0 and num2 == 0 and num3 == 0:
            break
        
        if num1 >= num2 and num1 >= num3:
            pos = [pos[0]-1, pos[1]-1]
            
            DNA_align = DNA_prot[pos[0]-1] + DNA_align
            prot_align = prot_seq[pos[1]-1] + prot_align
            
            max_num = num1
            
        elif num2 >= num3:
            pos = [pos[0], pos[1]-1]
            
            DNA_align = "-" + DNA_align
            prot_align = prot_seq[pos[1]-1] + prot_align
            
            max_num = num2
        
        else:
            pos = [pos[0]-1, pos[1]]
            
            DNA_align = DNA_prot[pos[0]-1] + DNA_align
            prot_align = "-" + prot_align
            
            max_num = num3
    
    align_DNA_all.append(DNA_align) #Add alignement to array
    align_prot_all.append(prot_align) #Add alignment to array

##############################################################################
# OUTPUT SCORE MATRIX AND SEQUENCE ALIGNMENT
##############################################################################

file = open("SubOpt_ScoreMatrix.csv", "w")

file.write(", ,")
for x in range(0, len(prot_seq)):
    file.write(prot_seq[x] + ",")
file.write("\n")

DNA_prot = " " + DNA_prot

for i in range(0, len(copy_matrix)):
    file.write((DNA_prot[i]) + ",")
    for j in range(0, len(copy_matrix[0])):
        file.write(str(copy_matrix[i][j]) + ",")
    file.write("\n")
file.close()

file2 = open("SubOpt_SequenceInfo.txt", "w")

file2.write("DNA:\t\t" + DNA_name[1:len(DNA_name)])
file2.write("Protein:\t" + prot_name[1:len(prot_name)])
file2.write("\nTranslated DNA: " + DNA_prot)
file2.write("\nProtein: " + prot_seq)
file2.write("\nUsed BLOSUM62 for substitution matrix\n")
file2.write("Gap Opening = " + str(gap_o))
file2.write("\nGap Extension = " + str(gap_e))

#Iterate through all alignments to print alignments between DNA and protein
for w in range(0, num_align):
    DNA_align_temp = align_DNA_all[w].replace("-", "")
    prot_align_temp = align_prot_all[w].replace("-", "")
    DNA_start = DNA_prot.index(DNA_align_temp)
    prot_start = prot_seq.index(prot_align_temp)
    if len(str(prot_start)) > len(str(DNA_start)):
        spaces = " " * (len(str(prot_start)) - len(str(DNA_start)))
        spaces2 = " " * (len(str(prot_start + len(prot_align_temp))) - len(str(DNA_start + len(DNA_align_temp))))
        align_DNA_all[w] = str(DNA_start) + spaces + "   " + align_DNA_all[w] + "   " + spaces2 + str(DNA_start + len(DNA_align_temp) - 1)
        align_prot_all[w] = str(prot_start + 1) + "   " + align_prot_all[w] + "   " + str(prot_start + len(prot_align_temp))
    else:
        spaces = " " * (len(str(DNA_start)) - len(str(prot_start)))
        spaces2 = " " * (len(str(DNA_start + len(DNA_align_temp))) - len(str(prot_start + len(prot_align_temp))))
        align_prot_all[w] = str(prot_start + 1) + spaces + "   " + align_prot_all[w] + "   " + spaces2 + str(prot_start + len(prot_align_temp))
        align_DNA_all[w] = str(DNA_start) + "   " + align_DNA_all[w] + "   " + str(DNA_start + len(DNA_align_temp) - 1)
    
    file2.write("\n\n\nALIGNMENT #" + str(w + 1) + "\n")
    file2.write("DNA:\t\t" + align_DNA_all[w] + "\n\t\t")
    for k in range(0, len(align_DNA_all[w])):
        if align_DNA_all[w][k] == "-" or align_prot_all[w][k] == "-" or align_DNA_all[w][k] != align_prot_all[w][k] or align_DNA_all[w][k] == " " or align_DNA_all[w][k].isnumeric() is True:
            file2.write(" ")
        else:
            file2.write("|")
    file2.write("\nProtein:\t" + align_prot_all[w] + "\n")


file2.write("\nNote: Look at SubOpt_ScoreMatrix.csv for full scoring matrix")

file2.close()



    
    
    
    
    
    