#reference values for calculations and conversions

#modules 


#dictionaries
	#3 letter and 1 letter codes (1 letter '2' 3 letter)
conversion = { 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','HYP':'O','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','GLP':'U','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
	#masses of AA residues
MASS = {'A': 71.08,'R': 156.19,'N': 114.11,'D': 115.09,'C': 103.15,'E': 129.12,'Q': 128.13,'G': 57.05,'H': 137.14,'O': 113.11,'I': 113.16,'L': 113.16,'K': 128.18,'M': 131.20,'F': 147.18,'P': 97.12,'U': 121.09,'S': 87.08,'T': 101.11,'W': 186.22,'Y': 163.18,'V': 99.13}

        #dissociation constant for -COOH group
PKA = {'A': 2.34,'R': 2.17,'N': 2.02,'D': 1.88,'C': 1.96,'E': 2.19,'Q': 2.17,'G': 2.34,'H': 1.82,'O': 1.82,'I': 2.36,'L': 2.36,'K': 2.18,'M': 2.28,'F': 1.83,'P': 1.99,'S': 2.21,'T': 2.09,'W': 2.83,'Y': 2.20,'V': 2.32}

        #dissociation constant for -NH3 group
PKB = {'A': 9.69,'R': 9.04,'N': 8.80,'D': 9.60,'C': 10.28,'E': 9.67,'Q': 9.13,'G': 9.60,'H': 9.17,'O': 9.65,'I': 9.60,'L': 9.60,'K': 8.95,'M': 9.21,'F': 9.13,'P': 10.60,'S': 9.15,'T': 9.10,'W': 9.39,'Y': 9.11,'V': 9.62}

        #dissociation constant for radical group
PKR = {'R': 12.48,'D': 3.65,'C': 8.18,'E': 4.25,'H': 6.00,'K': 10.53,'Y': 10.07}


#dictionary value test function (for comparing)
def referencecheck(dictionary):
    for a, b in dictionary.items():
        a = str(a)
        b = str(b)
        print(a + ' = ' + b)

