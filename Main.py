#external modules
import numpy as np
import matplotlib.pyplot as plt


#personal modules
    #reference values for calculations and conversions
import reference as ref
from reference import conversion
from reference import mass
from reference import pka
from reference import pkb
from reference import pkr

    #properties of chemical species
import properties as prop
from properties import positive
from properties import negative

#variables
    #x-axis values for pH plot, (decrease last number for higher resolution)
x = np.arange(0,14.25,0.25)

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
        temp = 0 #introduce temporary variable to add up equation for charges
        for i in range(int(len(self.sequence))): #calculate the charge on the radicals
            if self.sequence[i] in prop.positive:
                temp = temp + 1/(1+10**(x-ref.pkr[self.sequence[i]])) #equation for charge if AA is positive
                pass
            elif self.sequence[i] in prop.negative:
                temp = temp - (1*10**(x-ref.pkr[self.sequence[i]]))/(1+10**(x-ref.pkr[self.sequence[i]])) #equation for charge if AA is negative
                pass
            else:
                pass

        temp = temp + 1/(1+10**(x-ref.pkb[self.sequence[0]])) #taking into account the dissociation of the -NH3 group
        temp = temp - (1*10**(x-ref.pka[self.sequence[-1]]))/(1+10**(x-ref.pka[self.sequence[-1]])) #taking into account the dissociation of the -COOH group
        
        y = temp #assign temporary value to a permanent one
        print(y)
        return plt.plot(x,y, label=self.name) #returns y values to be plotted

        

#functions (mainly for interfacing)
    #determines how many sequences need to be analyzed
def inputnumber(): 
    global sequences
    number = int(input('How many sequences would you like to input?\n>>>'))
    sequences = list()
    for i in range(0,number):
        sequences.append(protein(input('name of protein ' + str(i+1) + ':'),input('sequence of protein ' + str(i+1) + ':')))
        sequences[i].chargeplot()

#interface
plt.figure(1)

print(x)

inputnumber()

#plt.plot(x,test2.chargeplot(),'b-', label=input('second inputted protein:'))
plt.xlabel('pH')
plt.ylabel('relative charge')
plt.legend() #get legend to show on plot

    #configure axes
ax = plt.gca() #get current axes (then spines to configure axes)
ax.spines['top'].set_color('none') #removes top black line on plot
ax.spines['bottom'].set_position('zero') #ensures x axis always crosses at y = 0
ax.spines['right'].set_color('none') #removes right black line on plot

    #x axis ticks  
start,end = ax.get_xlim() #gets min and max of x axis
end = end+1
ax.set_xticks(np.arange(0,round(end),1)) #sets major ticks on x axis
ax.set_xticks(np.arange(round(start),round(end),0.25),True) #sets minor ticks on x axis

    #y axis ticks 
start,end = ax.get_ylim() #gets min and max of y axis
end = end + 1
#ax.set_yticks(np.arange(round(start),round(end),5)) #sets major ticks on y axis
ax.set_yticks(np.arange(round(start),round(end),1),True) #sets minor ticks on y axis

plt.grid()

plt.show()




#NSDSECPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWELR

#MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN

#MAAGSITTLPALPEDGGSGAFPPGHFKDPKRLYCKNGGFFLRIHPDGRVDGVREKSDPHIKLQLQAEERGVVSIKGVCANRYLAMKEDGRLLASKCVTDECFFFERLESNNYNTYRSRKYTSWYVALKRTGQYKLGSKTGPGQKAILFLPMSAKS