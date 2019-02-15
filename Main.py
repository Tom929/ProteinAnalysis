#external modules
import numpy as np
import matplotlib.pyplot as plt


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
        self.mass = int(temp)/1000 #conversion to kDA
        return plt.bar(self.name,self.mass,label=str(self.mass) + ' kDA') #returns plot of kDA

#functions (mainly for interfacing)
    #determines how many sequences need to be analyzed
def inputnumber(): 
    global sequences
    sequences = list()
    number = int(input('How many sequences would you like to input?\n>>>'))
    for i in range(0,number):
        #sequences.append(protein(input('name of protein ' + str(i+1) + ':'),input('sequence of protein ' + str(i+1) + ':'))) #adds proteins to sequences list
        sequences.append(data.access_sequence(input('Protein ' + str(i+1)+ ': '))) #list of sequences accessed at SwissProt
    return

def generatechargeplot():
    chargeplot = plt.figure()
    x = np.arange(0,14.25,0.25) 
    for i in range(0,int(len(sequences))): #generates the y values
        sequences[i].chargeplot() 

    plt.grid() # adds gridlines

    plt.xlabel('pH') 
    plt.ylabel('Relative Charge')
    plt.legend() 

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
    ax.set_yticks(np.arange(round(start),round(end),1),True) #sets minor ticks on y axis
    
    return chargeplot 

def generatemassplot():
    massplot = plt.figure()
    for i in range(0,int(len(sequences))): #generates the y values
        sequences[i].massplot()

    plt.xlabel('Protein')
    plt.ylabel('Mass (kDa)')
    plt.legend()

    #configure axes
    bx = plt.gca() #get current axes (then spines to configure axes)
    bx.spines['top'].set_color('none') #removes top black line on plot
    bx.spines['bottom'].set_position('zero') #ensures x axis always crosses at y = 0
    bx.spines['right'].set_color('none') #removes right black line on plot
    bx.set_axisbelow(True)

    #y axis ticks 
    start,end = bx.get_ylim() #gets min and max of y axis
    start = 0
    end = end + 1
    bx.set_yticks(np.arange(0,round(end),5)) #sets major ticks on x axis
    bx.set_yticks(np.arange(round(start),round(end),1),True) #sets minor ticks on y axis
    bx.yaxis.grid(which='major')

    return massplot

def filecreate(name):
    file = open(name,'w')

    for i in range(len(sequences)):
        file.write(sequences[i].name + ':\n' +sequences[i].sequence + '\n\n')
    
    file.close()

    return



#interface
inputnumber() # user input of protein sequences

print(sequences)

#filecreate('current')

#generatechargeplot()
#generatemassplot()



#plt.show()





#Some sequences to test:

    #EGF
#NSDSECPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWELR

    #Insulin
#MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN

    #FGF2
#MAAGSITTLPALPEDGGSGAFPPGHFKDPKRLYCKNGGFFLRIHPDGRVDGVREKSDPHIKLQLQAEERGVVSIKGVCANRYLAMKEDGRLLASKCVTDECFFFERLESNNYNTYRSRKYTSWYVALKRTGQYKLGSKTGPGQKAILFLPMSAKS

    #TGFB1
#MPPSGLRLLLLLLPLLWLLVLTPGRPAAGLSTCKTIDMELVKRKRIEAIRGQILSKLRLASPPSQGEVPPGPLPEAVLALYNSTRDRVAGESAEPEPEPEADYYAKEVTRVLMVETHN
