import numpy as np

#define matrices for standard gates
X = np.array([[0.0,1.0],[1.0,0.0]])
Y = np.array([[0,-1.0j],[1.0j,0.0]])
Z = np.array([[1.0,0.0],[0.0,-1.0]])
I = np.eye(2)
H = np.array([[1.0,1.0],[1.0,-1.0]])*(1/np.sqrt(2.0))
CNOT = np.array([[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,0.0,1.0],[0.0,0.0,1.0,0.0]])
CZ = np.array([[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,-1.0]])

gate_dict = {
    "X": X, 
    "Y": Y,
    "Z": Z, 
    "I": I, 
    "H": H, 
    "CNOT": CNOT, 
    "CZ": CZ
}

class Gate: 
    def __init__(self, name, qubits, angles=[]):
        self.name = name
        self. angles = angles
        self.qubits = qubits

    def matrix(self):
        if self.name in gate_dict:
            return gate_dict[self.name]
        elif (self.name == 'RZ'):
            return np.array([[np.cos(self.angles[0]/2)-1.0j*np.sin(self.angles[0]/2), 0],[0, np.cos(self.angles[0]/2)+1.0j*np.sin(self.angles[0]/2)]])
        elif (self.name == 'RX'):
            return np.array([[np.cos(self.angles[0]/2), -1.0j*np.sin(self.angles[0]/2)],[-1.0j*np.sin(self.angles[0]/2), np.cos(self.angles[0]/2)]])
        elif (self.name == 'RY'):
            return np.array([[np.cos(self.angles[0]/2), -np.sin(self.angles[0]/2)],[np.sin(self.angles[0]/2), np.cos(self.angles[0]/2)]])
        else:
            print("Error: ", self.name, " is not a known gate name!")
            exit()

    def gate_from_Pauli(self, pauli):
        self.name = pauli.name
        self.qubits = pauli.qubit
        self.angles = []

    def print_gate(self):
        if (self.angles != []):
            print(self.name,"(",self.angles[0],")",self.qubits) 
        else: 
            print(self.name,self.qubits)  

class Pauli: 
    def __init__(self, name, qubit):
        self.name = name
        self.qubit = qubit

class Term:
    def __init__(self, paulis, coeff):
        self.paulis = paulis
        self.coeff = coeff

class Hamiltonian:
    def __init__(self, nqubits, terms):
        self.nqubits = nqubits
        self.terms = terms

    def add_term(self, term):
        self.term.append(term)

    def matrix(self):
        dim = 2**self.nqubits
        ham_mat = np.zeros((dim,dim))
        for term in self.terms:
            kron_list = ['I'] * self.nqubits
            for pauli in term.paulis:
               kron_list[pauli.qubit] = pauli.name
            for q in range(self.nqubits-1):
                if (q==0):
                    mat = np.kron(gate_dict[kron_list[0]], gate_dict[kron_list[1]])
                else:          
                    mat = np.kron(mat, gate_dict[kron_list[q+1]])
            mat = term.coeff * mat 
            ham_mat = ham_mat + mat
        return ham_mat
    
    def show(self):
        dim = 2**self.nqubits
        ham_mat = self.matrix()
        for i in range(dim):
            for j in range(dim):
                print(ham_mat[i][j])



class Program:
    def __init__(self, nqubits):
        self.nqubits = nqubits
        self.gates = [] 

    def add_instr(self, gate_list):
        for gate in gate_list:
            self.gates.append(gate)

    def get_U(self):
        dim = 2**self.nqubits
        matU = np.eye(dim)
        for gate in self.gates:
            #make sure gate matrix has dimension of system
            if(self.nqubits > 1):
                kron_list = [I] * self.nqubits
                kron_list[gate.qubits[0]] = gate.matrix()
                num_prods = self.nqubits-1
                if (len(gate.qubits) == 2):
                    num_prods = self.nqubits-2
                    if (self.nqubits > 2):
                        #for a 2-qubit gate, the last identity should be removed from kron_list
                        kron_list.pop()
                    else: mat = gate.matrix()
                for q in range(num_prods):
                    if (q==0):
                        mat = np.kron(kron_list[0], kron_list[1])
                    else:
                        mat = np.kron(mat, kron_list[q+1])
            else: mat = gate_dict[gate.matrix]
            matU = np.matmul(matU,mat)
        return matU  
