import sys
import pandas as pd

# creating dictionary of monoisotopic averages from the Mascot database search Amino acid reference data (2021)
Mmass = {
    "A" : 71.0371, "C" : 103.0092, "D" : 115.0269, "E" : 129.0426,
    "F" : 147.0684, "G" : 57.0215, "H" : 137.0589, "I" : 113.0841,
    "K" : 128.0950, "L" : 113.0841, "M" : 131.0405, "N" : 114.0429,
    "P" : 97.0528, "Q" : 128.0586, "R" : 156.1011, "S" : 87.0320,
    "T" : 101.0477, "V" : 99.0684, "W" : 186.0793, "Y" : 163.0633,
    "\s" : 0.0, "*" : 0.0
}

#Creating a calculator for monoisotopic mass values 
MH2O = 18.0106            # the mass of the terminating groups for monoisotopic masses
def Mcalculator(sequence): 
    values = []           #creating a list of average mass values
    for character in sequence:
        for key, value in Mmass.items(): #returns a view object containing the key-value pairs of the dictionary, as tuples in a list.
            if key == character:         #if the key is a character
                values.append(value)     #add the values of the keys to the list values
    return("{:.4f}".format(sum(values) + MH2O)) #return the sum of the values in the list along with the mass of water to 4 decimal places

# creating a dictionary of average masses from the Mascot database search Amino acid reference data (2021)
Amass = {
    "A" : 71.0779, "C" : 103.1429, "D" : 115.0874, "E" : 129.114,
    "F" : 147.1739, "G" : 57.0513, "H" : 137.1393, "I" :113.1576,
    "K" : 128.1723, "L" : 113.1576, "M" : 131.1961, "N":114.1026,
    "P" : 97.1152, "Q" : 128.1292, "R" : 156.1857, "S" : 87.0773,
    "T" : 101.1039, "V" : 99.1311, "W" : 186.2099, "Y" :163.1733,
    "\s" : 0.0, "*" : 0.0
}

#creating a calculator for the average mass values 
AH2O =  18.0153           #the mass of the terminating groups for average masses 
def Acalculator(sequence):
    values = []           #creating a list of average mass values
    for character in sequence:
        for key, value in Amass.items(): #returns a view object containing the key-value pairs of the dictionary, as tuples in a list.
            if key == character:         #if the key is a character
                values.append(value)     #add the values of the keys to the list values
    return("{:.4f}".format(sum(values) + AH2O)) #return the sum of the values in the list along with the mass of water to 4 decimal places 

# creating a dicitionary of the hydrophobicity index (Argos et al., 1982)
HI = {
    "A" : 0.61, "C" : 1.07, "D" : 0.46, "E" : 0.47,
    "F" : 2.02, "G" : 0.07, "H" : 0.61, "I" : 2.22,
    "K" : 1.15, "L" : 1.53, "M" : 1.18, "N": 0.06,
    "P" : 0.05, "Q" : 0.00, "R" : 0.60, "S" : 0.05,
    "T" : 0.05, "V" : 1.32, "W" : 2.65, "Y" : 1.88,
    "\s" : 0.0, "*" : 0.0
}

# creating a calculator for the hydrophobicity index
def HIcalculator(sequence):
    values = []               #creating a list of hydrophobicity values
    for character in sequence:
        for key, value in HI.items(): #returns a view object containing the key-value pairs of the dictionary, as tuples in a list.
            if key == character:      #if the key is a character
                values.append(value)  #add the values of the keys to the list values
    return("{:.4f}".format((sum(values))/(len(values)))) #return the sum of the values in the list to 4 decimal places


#creating a list of names 
def listNames(fileName):
    file = open(fileName, 'r')
    seqNames = []  #  a list containing all the sequence names 
    
    
    for line in file:
        if line.startswith(" >"):
            data = line.split()  #splitting this line into a list
            name_with_arrow = data[0] #extract 1st element from the list 
            name = name_with_arrow [1:] #removing the arrow 
            if name:
                seqNames.append(name) # spill over section - add last seq onto list 
            else:
                seqNames[name] += line.rstrip('\n') #remove the new trailing line 
    file.close()
    return (seqNames)


file = sys.argv[1] #the file is put in on the command line - will be one of the four files generated from task 2
names = listNames(file)

# create list of peptides
def listPeptides(fileName):
    file = open(fileName, 'r')
    peptides = []  #  a list containing the peptide numbers 
    
    
    for line in file:
        if line.startswith(" >"):
            data = line.split() #splitting first line from string into a list
            peptide = data[2]   #extract the 3rd element from the list
            if peptide:
                peptides.append(peptide) # spill over section add last peptide number onto list 
            else:
                peptides[peptide] += line.rstrip('\n') #remove the new trailing line
    file.close()
    return (peptides)


file = sys.argv[1]
peptide = listPeptides(file)

# creating a list of the amino acid sequences 
def listSeqs (fileName):
    fileObj = open(fileName, 'r')
    sequences = []   #  a list, to contain all our sequences
    seqFrags = []    #  a list, to contain local lines of bits of a sequence
    seq = ''         #  initialise our sequence placeholder to blank 
    
    for line in fileObj:
        if line.startswith(' >'):     # a new sequence?
            if seqFrags:             # Have we already read seq data in?
                seq = ''.join(seqFrags)    # join is more efficient
                sequences.append(seq)
            seqFrags = []           # reset our fragments list
        else:
            seq = line.rstrip()  # remove newline, concat to growing seq string
            seqFrags.append(seq)
    if seqFrags:
        seq = ''.join(seqFrags)
        sequences.append(seq) # spill over section - add last seq onto list    
    fileObj.close()
    return sequences


file = sys.argv[1]
mySeqs = listSeqs(file)


p = ['1'] * len(mySeqs) #adding 1 to all peptide ion values

#creating a list of the digestive enzyme used to create the peptide sequences corresponding to the input FASTA file generated from task 2 
if sys.argv[1] == 'trypsin (1).fasta': 
    enzyme = ['Trypsin'] * len(mySeqs) 
if sys.argv[1] == 'lys-c (1).fasta':
    enzyme = ['Endoproteinase_Lys-C'] * len(mySeqs)
if sys.argv[1] == 'glu-c (1).fasta':
    enzyme = ['V8_proteinase_(Glu-C)'] * len(mySeqs)
if sys.argv[1] == 'arg-c (1).fasta':
    enzyme = ['Endoproteinase_Arg-C'] * len(mySeqs)

#creating a list of all monoisotopic masses 
Mmasses = [] 
for i in range(0,len(mySeqs)):
    Mmasses.append(Mcalculator(mySeqs[i])) #add the monoisotopic masses of the sequences in mySeqs into list 

#creating a list of all the average masses 
Amasses = []  
for i in range(0,len(mySeqs)):
    Amasses.append(Acalculator(mySeqs[i])) #add the average masses of the sequences in mySeqs into list 

# creating a list of all the hydrophobicites 
Hydro = [] 
for i in range(0,len(mySeqs)):
    Hydro.append(HIcalculator(mySeqs[i])) #add the hydrophobicities of the sequences in mySeqs into list 



#writing a new file to store the pep-masses 
infile = sys.argv[1]
outfile = 'pep-masses.txt'

#generating pep-mass files for each enzyme so that the output file will correspond to the input enzyme FASTA file and the output file will not be overwritten 
if sys.argv[1] == 'trypsin (1).fasta':         #if 'trypsin (1).fasta' written on the command line 
    fh = open('pep-masses(trypsin).txt', 'w')  #the data generated will be written to a txt file named 'pep-masses(trypsin).txt'
if sys.argv[1] == 'lys-c (1).fasta':
    fh = open('pep-masses(lys).txt', 'w')
if sys.argv[1] == 'glu-c (1).fasta':
    fh = open('pep-masses(glu).txt', 'w')
if sys.argv[1] == 'arg-c (1).fasta':
    fh = open('pep-masses(arg).txt', 'w')

#giving the user the choice whether to view the monoisotopic masses or average masses 
mass = input('what type of mass-to-charge ratios would you like to view? monoisotopic or average?')

#creating a table with all the data in using pandas 
#writing the table with the monoisotopic masses to the ouput file 
if mass =="monoisotopic":
    print(pd.DataFrame(list(zip(names, enzyme, peptide, p, Mmasses, Hydro, mySeqs)),
            columns=['Name', 'Enzyme', 'Peptide', 'p', 'mass-to-charge', 'Hydrophobicity Index','Sequence']), pd.set_option('max_rows', None), pd.reset_option("expand_frame_repr"), file = fh)
#pd set option allows all the rows to be printed to the file 

#writing the table with the average masses to the output file 
if mass =="average":
    print(pd.DataFrame(list(zip(names, enzyme, peptide, p, Amasses, Hydro, mySeqs)),
            columns=['Name', 'Enzyme', 'Peptide', 'p', 'mass-to-charge', 'Hydrophobicity Index','Sequence']), pd.set_option('max_rows', None), pd.reset_option("expand_frame_repr"), file = fh)