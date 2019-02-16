#external modules
import numpy as np
import matplotlib.pyplot as plt
from Bio import ExPASy
from Bio import SwissProt


#personal modules
    #reference values for calculations and conversions
import reference as ref
from reference import conversion
from reference import MASS
from reference import PKA
from reference import PKB
from reference import PKR

    #properties of chemical species
import properties as prop
from properties import POSITIVE
from properties import NEGATIVE

    #database access
import database as data

   #figure plotting
import graph

#classes (where the magic happens)

    #attributes of the inputted polypeptide seqeuence
class protein:
    def __init__(self,name,sequence):
        self.name = name
        self.sequence = sequence.upper() #sequence is capitalized to standardize inputs
        return

        #converts three letter code to one letter code
    def convert(self): 
        temp = '' #introduce temporary variable to concatenate full one letter code
        for i in range(1,int(len(self.sequence)/3+1)):
            #print(self.sequence[3*(i-1):3*(i-1)+3]) #only include when debugging!
            temp = temp + conversion[self.sequence[3*(i-1):3*(i-1)+3]]
        self.sequence = temp
        return

        #calculates relative charge on peptide at specific pH
    def chargeplot(self):
        x = np.arange(0,14.25,0.25) #x values for plot and calculation of charge
        temp = 0 #introduce temporary variable to add up equation for charges
        for i in range(int(len(self.sequence))): #calculate the charge on the radicals
            if self.sequence[i] in prop.POSITIVE:
                temp = temp + 1/(1+10**(x-ref.PKR[self.sequence[i]])) #equation for charge if AA is positive
                pass
            elif self.sequence[i] in prop.NEGATIVE:
                temp = temp - (1*10**(x-ref.PKR[self.sequence[i]]))/(1+10**(x-ref.PKR[self.sequence[i]])) #equation for charge if AA is negative
                pass
            else:
                pass

        temp = temp + 1/(1+10**(x-ref.PKB[self.sequence[0]])) #taking into account the dissociation of the -NH3 group
        temp = temp - (1*10**(x-ref.PKA[self.sequence[-1]]))/(1+10**(x-ref.PKA[self.sequence[-1]])) #taking into account the dissociation of the -COOH group
        
        self.charge = temp #assign new attribute for charge at different pH
        return plt.plot(x,self.charge, label=self.name) #returns plot of charge
    
        #calculates the approximate masses of proteins in kDa (at physiological pH)
    def massplot(self):
        temp = 0
        for i in range(int(len(self.sequence))): #add up the masses of the residues in Da
            temp = temp + ref.MASS[self.sequence[i]]
        temp = temp + 16 #include the extra oxygen present in the last carboxy group
        self.mass = int(temp)/1000 #conversion to kDA (conversion to int gives 3 d.p.)
        return plt.bar(self.name,self.mass,label=str(self.mass) + ' kDA') #returns plot of kDA

#functions (mainly for interfacing)
    #determines how many sequences need to be analyzed
def inputnumber(input): 
    global sequences
    sequences = list()
    number = int(input('How many sequences would you like to input?\n>>>'))
    for i in range(0,number):
        temp = input('Protein ' + str(i+1)+ ': ')
        sequences.append(protein(data.access_name(temp), data.access_sequence(temp))) #list of sequences accessed at SwissProt
        print(data.access_sequence(temp))
    return



def filecreate(name):
    file = open(name,'w')

    for i in range(len(sequences)):
        file.write(sequences[i].name + ':\n' +sequences[i].sequence + '\n\n')
    
    file.close()

    return



#interface
#inputnumber() # user input of protein sequences

#print(len(sequences))
sequences = list()
file = open('ProteinInput.txt','r')

input = file.read().splitlines()

for i in input:
    sequences.append(protein(data.access_name(i), data.access_sequence(i))) #list of sequences accessed at SwissProt
    print(data.access_name(i))
    print(data.access_sequence(i))

file.close()



#filecreate('current') ## in progress to write text file with current inputted sequences

graph.generatechargeplot(sequences)
graph.generatemassplot(sequences)

plt.show()






