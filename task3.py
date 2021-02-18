import numpy as np
import random
import json
from math import log2


#Define the ground state for n_quibits:

def get_ground_state(n_qubits):
    qc=[1]
    for i in range(0,2**(n_qubits)-1):
        qc.append(0);
    return  qc


#Define 1-quibit Gates

#X-Not Gate
X=np.array([
[0,1],
[1,0]
])

#Hadamard gate
H = np.array([
[1/np.sqrt(2), 1/np.sqrt(2)],
[1/np.sqrt(2), -1/np.sqrt(2)]
])

#Z gate
Z=np.array([
[1,0],
[0,-1]
])

#Define function to get the operators to apply the gates in the circuit 

def get_operator(n_qubits,unitary, target_qubits):

    #Define 2x2 Identity
    I=np.identity(2)

    #Evaluate which gate have been selected to apply

    #Two qubits Control-X gate

    if (unitary=='CX'):

        #Define Pojectors
        #|0><0|
        P0x0=np.array([
        [1,0],
        [0,0]
        ])
        #|1><1|
        P1x1=np.array([
        [0,0],
        [0,1]
        ])

        #Define the tensor products order here I define two producs |0><0|(x)I and |1><1|(x)X
        #Need to apply the projectors and the gate in order to get the desired Control-X gate

        #Using auxiliar matrices E11 E12 E21 E22 to storage the partical values of |0><0|(x)I and |1><1|(x)X operations

        #Determine if first qubit is the control one on first term |0><0|(x)I

        if (target_qubits[0]==0):

            E11=P0x0     #If first qubit is the control operation starts with P0x0

        else:

            E11=I        #Othewise operations starts with 2x2 Identity matix
            
        #Determine if first qubit is the control or target one on second term |1><1|(x)I

        if (target_qubits[0]==0):

            E12=P1x1     #If first qubit is the control operation starts with P1x1

        elif(target_qubits[1]==0):

            E12=X        #If fist qubit is the target operation starts with the X-gate

        else:

            E12=I        #Otherwise operarion starts with 2x2 Identity matrix

        #Continue Looking for control and target qubits for the n-1 last qubits, and do the Kronecker products to get the final nxn matrix operator

        for i in range(1,n_qubits):

            if (i==target_qubits[0]): #Looking for control qubit on first term |0><0|(x)I 
                E21=P0x0
            else:
                E21=I

            if (i==target_qubits[0]): #Looking for control qubit on second term |1><1|(x)I
                E22=P1x1
            elif (i==target_qubits[1]):
                E22=X                 #Loking for target qubit on second term |1><1|(x)I
            else:
                E22=I

            #Apply the Kronecker product to single qubit gates 
            E11=np.kron(E11,E21)
            E12=np.kron(E12,E22)

        #get the Controlled-n operator 
        O=E11+E12
    else:

        #For 1-quibit gates

        #Possible 1-qubit gates to apply
        my_gate= {
                "X":X,
                "H":H,
                "Z":Z
        }

        #Define the 1-qubit gate operatos doing Kronecker product with 2x2 Identity matrices and the 2x2 matrix that represents the gate

        #Find target qubit in which will apply the unitary gate

        #Using auxilar variables P1 and P2 to storage partial values of the operator

        if (target_qubits[0]==0):
            P1=my_gate.get(unitary) #If targer qubit is the first one use the gate's matrix representation
        else:
            P1=I                    #Otherwise use Identity matrix

        #Look for target qubit on subsequent qubits

        for i in range(1,n_qubits):
            if (i==target_qubits[0]):
                P2=my_gate.get(unitary) #Looking for target qubit 
            else:
                P2=I

            #Apply the Kronecker product to sibgle qubit gates
            O=np.kron(P1,P2)

            #Replace P1 value to continue the operation until we get the desired matrix representation of the gate for a n_qubits state
            P1=O
    return O

#Function to apply a gate on the initial ground state and get the final state
def run_program(initial_state,program):

    #get qubits number
    nq=int(log2(len(initial_state)))

    #Create a list to storage all the gates selected to apply
    gates=[] 

    #Use get_operator function to define the matrices to apply the gates on intial state
    for op in program:
        gates.append(get_operator(nq,op['gate'],op['target']))         

   #given thatt the gates are determined in left to right order, it is necessary  to reverse the gates list to apply them on initial state in correct order which is right to left
    gates=gates[::-1]
    
    #Apply the Gates on initial_state to get final_state

    for i in range(len(gates)): #Apply dot product with a selected gate and the initial state to get a new state
        
        final_state=np.dot(gates [i],initial_state)
        initial_state=final_state

    return final_state 

#This funtion produces a list of all possibles binary states
def state_list(n_qubits):
    states=[]
    for i in range(2**n_qubits):
        s=bin(i)[2:] #binary convertion of i
        states.append(str(s.zfill(n_qubits)))
    return states

#Funtion to make a measurement on the quantum circuit

def measure_all(states,state_vector):

    prob=np.abs(state_vector)**2 #given a state_vector with probabily ampitudes, determine the probabilities to get any state after performing a measure on the circuit

    measure=random.choices(states,weights=prob,k=1) #weighed random selection on the possibles outcomes from the circuit

    ms=str(measure)[2:-2] #String giving the outcome state

    j=states.index(ms) #The variable j the index of the state on states list

    return j

#Do num_shots measures

def get_counts(state_vector, num_shots):

    nq=int(log2(len(state_vector))) #obtain the number of qubits

    times=[] #list to storage the total times that every possible outcome is obtained 

    for i in range(len(state_vector)):
        times.append(0)


    states=state_list(nq) #Create a list of all possible outcome states 

    #Perform num_shots measures on the circuit
    for i in range(num_shots):

        j=measure_all(states,state_vector)

        times[j]+=1 #add +1 to the number of appearances of state j

    #filter the states to return only those obtained at least one time    
        fs=[]
        ft=[]

    for i in range(len(times)):
        if (times[i]!=0):
            fs.append(states[i])
            ft.append(times[i]) 


    #Create a Dictornary with all the states obtained at least once and its total outcome times
    counts=dict(zip(fs,ft))
    return counts


#Test Program

#Initialize ground_state
my_qpu=get_ground_state(8)

#Configure program selecting the gates to apply
my_cir=[{"gate":"H","target":[0]}]

#num_shots
num_shots=1000

#Apply gates to get final state
f_state=run_program(my_qpu,my_cir)

#Measure for num_shots counts
counts=get_counts(f_state,num_shots)

#Print the results
print(json.dumps(counts,indent=2))





