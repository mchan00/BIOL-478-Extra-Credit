# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:17:09 2021

@author: Matthew Chan
"""
'''
Look at LocalAlign_code for more comments since it it quite similar
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

#First row can go negative by adding gap penalties
first_row = [0]
for x in range(0, len(prot_seq)):
    first_row.append(x * gap_e + gap_o)
score_matrix.append(first_row)

#Determines gap penalty depending on the column number
for y in range(0, len(DNA_prot)):
    temp = [(y)*gap_e + gap_o]
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
        
        score_matrix[x+1][y+1] = max_val

##############################################################################
# TRACEBACK
##############################################################################

pos = [len(DNA_prot),len(prot_seq)]

DNA_align = ""
prot_align = ""

#max_num = max_index[2]

while pos[0] != 0 and pos[1] != 0:
    
    num1 = score_matrix[pos[0]-1][pos[1]-1]
    num2 = score_matrix[pos[0]][pos[1]-1]   #up?
    num3 = score_matrix[pos[0]-1][pos[1]]   #left?

    score_matrix[pos[0]][pos[1]] = "*" + str(score_matrix[pos[0]][pos[1]])
    
    if (pos[0] - 1) == 0 and (pos[1] - 1) == 0:
        break
    
    #End traceback is different than rest of the traceback
    if pos[0] == len(DNA_prot) or pos[1] == len(prot_seq):
        if num1 >= num2 and num1 >= num3:
            
            DNA_align = DNA_prot[pos[0]-1] + DNA_align
            prot_align = prot_seq[pos[1]-1] + prot_align
            pos = [pos[0]-1, pos[1]-1]
            
        elif num2 >= num3:
            
            DNA_align = "-" + DNA_align
            prot_align = prot_seq[pos[1]-1] + prot_align
            pos = [pos[0], pos[1]-1]
            
        else:
            
            DNA_align = DNA_prot[pos[0]-1] + DNA_align
            prot_align = "-" + prot_align
            pos = [pos[0]-1, pos[1]]
    
    elif num1 >= num2 and num1 >= num3:
        
        pos = [pos[0]-1, pos[1]-1]
        
        DNA_align = DNA_prot[pos[0]-1] + DNA_align
        prot_align = prot_seq[pos[1]-1] + prot_align
        
        
    elif num2 >= num3:
        
        pos = [pos[0], pos[1]-1]
        
        DNA_align = "-" + DNA_align
        prot_align = prot_seq[pos[1]-1] + prot_align
        
    else:
        
        pos = [pos[0]-1, pos[1]]
        
        DNA_align = DNA_prot[pos[0]-1] + DNA_align
        prot_align = "-" + prot_align
        

##############################################################################
# OUTPUT SCORE MATRIX AND SEQUENCE ALIGNMENT
##############################################################################

file = open("Global_ScoreMatrix.csv", "w")

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

file2 = open("Global_SequenceInfo.txt", "w")

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
file2.write("\nNote: Look at Global_ScoreMatrix.csv for full scoring matrix")

file2.close()





    
    
    
    
    
    