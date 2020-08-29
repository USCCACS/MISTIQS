#import necessary libraries
import numpy as np
from quantum_circuits import Program, Gate
from ds_compiler import ds_compile
import os




#Create Data directory
current=os.getcwd()
newdir="Data"
path = os.path.join(current, newdir) 
if not os.path.isdir(path):
    os.makedirs(path)



class Heisenberg:

    def __init__(self,file="input_file.txt",log="logfile.txt"):
        #%matplotlib inline 
        input_file=open(file,'r')
        data=input_file.readlines()
        completename = os.path.join(path,log)
        self.namevar=str(completename)
        with open(self.namevar,'w') as tempfile:
            tempfile.write("***MISTIQS Session Log File***\n\n")
        self.H_BAR = 0.658212    # eV*fs

        #Default Parameters
        self.JX=self.JY=self.JZ=self.h_ext=0
        self.ext_dir="Z"
        self.num_qubits=2
        self.initial_spins="1,1"
        self.delta_t=1
        self.steps=1
        self.QCQS="QS"
        self.shots=1024
        self.noise_choice="n"
        self.device_choice="ibmq_rome"
        self.plot_flag="y"
        self.freq=0
        self.time_dep_flag="n"
        self.custom_time_dep="n"
        self.print_bool=0 #controls print statements for domain specific compiler integration
        self.smart_bool=False #controls transpiling in presence of domain specific compilation
        self.circuits_list=[]
        self.backend="ibm"
        self.ibm_circuits_list=[]
        self.rigetti_circuits_list=[]
        self.cirq_circuits_list=[]
        self.auto_smart_compile="y"
        self.default_compiler="native" #native or domain specific
        self.compile="y"

        from numpy import cos as cos_func
        self.time_func=cos_func

        for i in range(len(data)-1):
            value=data[i+1].strip()
            if "*JX" in data[i]:
                self.JX=float(value)
            elif "*JY" in data[i]:
                self.JY=float(value)
            elif "*JZ" in data[i]:
                self.JZ=float(value)
            elif "*h_ext" in data[i]:
                self.h_ext=float(value)
            elif "*initial_spins" in data[i]:
                self.initial_spins=value
            elif "*delta_t" in data[i]:
                self.delta_t=int(value)
            elif "*steps" in data[i]:
                self.steps=int(value)
            elif "*num_qubits" in data[i]:
                self.num_qubits=int(value)
            elif "*QCQS" in data[i]:
                self.QCQS=value
            elif "*device" in data[i]:
                self.device_choice=value
            elif "*backend" in data[i]:
                self.backend=value
            elif "*noise_choice" in data[i]:
                self.noise_choice=value
            elif "*plot_flag" in data[i]:
                self.plot_flag=value
            elif "*shots" in data[i]:
                self.shots=int(value)
            elif "*freq" in data[i]:
                self.freq=float(value)
            elif "*time_dep_flag" in data[i]:
                self.time_dep_flag=value
            elif "*default_compiler" in data[i]:
                self.default_compiler=value
            elif "*compile" in data[i]:
                self.compile=value
            elif "*ext_dir" in data[i]:
                self.ext_dir=value
            elif "*auto_smart_compile" in data[i]:
                self.auto_smart_compile=value
            elif "*custom_time_dep" in data[i]:
                self.custom_time_dep=value
                if self.custom_time_dep in "y":
                    from time_dependence import external_func
                    print("Found an external time dependence function")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Found an external time dependence function\n")
                    self.time_func=external_func


        self.initial_spins=self.initial_spins.split(',')



        self.total_time=int(self.delta_t*self.steps)



        self.flip_vec=np.zeros(self.num_qubits)
        index=0
        for spin in self.initial_spins:
            if int(spin)==-1:
                self.flip_vec[index]=1
                index+=1
            elif int(spin)==1:
                self.flip_vec[index]=0
                index+=1
            else: 
                print('Invalid spin entered')
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Invalid spin entered\n")


    def local_evolution_circuit(self,evol_time): #creates evolution circuit in local program
    #Initial flipped spins are not implemented in this function due to the need for "barrier". Need to do that outside of this.
        prop_steps = int(evol_time / self.delta_t)  # number of propagation steps
        P=Program(self.num_qubits)
        for step in range(prop_steps):
            t = (step + 0.5) * self.delta_t
            if "n" in self.time_dep_flag:
                psi_ext = -2.0 * self.h_ext *self.delta_t / self.H_BAR
            elif "y" in self.time_dep_flag:
                if "y" in self.custom_time_dep:
                    psi_ext = -2.0 * self.h_ext * self.time_func(t)*self.delta_t / self.H_BAR
                elif "n" in self.custom_time_dep:
                    psi_ext=-2.0*self.h_ext*np.cos(self.freq*t)*self.delta_t/self.H_BAR
                else:
                    print("Invalid selection for custom_time_dep parameter. Please enter y or n.")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Invalid selection for custom_time_dep parameter. Please enter y or n.\n")
                    break
            ext_instr_set=[]
            XX_instr_set=[]
            YY_instr_set=[]
            ZZ_instr_set=[]
            for q in range(self.num_qubits):
                if self.ext_dir in "X":
                    ext_instr_set.append(Gate('RX', [q], angles=[psi_ext]))
                elif self.ext_dir in "Y":
                    ext_instr_set.append(Gate('RY', [q], angles=[psi_ext]))
                elif self.ext_dir in "Z":
                    ext_instr_set.append(Gate('RZ', [q], angles=[psi_ext]))
            psiX=-2.0*(self.JX)*self.delta_t/self.H_BAR
            psiY=-2.0*(self.JY)*self.delta_t/self.H_BAR
            psiZ=-2.0*(self.JZ)*self.delta_t/self.H_BAR

            for q in range(self.num_qubits-1):
                XX_instr_set.append(Gate('H',[q]))
                XX_instr_set.append(Gate('H',[q+1]))
                XX_instr_set.append(Gate('CNOT',[q, q+1]))
                XX_instr_set.append(Gate('RZ', [q+1], angles=[psiX]))
                XX_instr_set.append(Gate('CNOT',[q, q+1]))
                XX_instr_set.append(Gate('H',[q]))
                XX_instr_set.append(Gate('H',[q+1]))

                YY_instr_set.append(Gate('RX',[q],angles=[-np.pi/2]))
                YY_instr_set.append(Gate('RX',[q+1],angles=[-np.pi/2]))
                YY_instr_set.append(Gate('CNOT',[q, q+1]))
                YY_instr_set.append(Gate('RZ', [q+1], angles=[psiY]))
                YY_instr_set.append(Gate('CNOT',[q, q+1]))
                YY_instr_set.append(Gate('RX',[q],angles=[np.pi/2]))
                YY_instr_set.append(Gate('RX',[q+1],angles=[np.pi/2]))

                ZZ_instr_set.append(Gate('CNOT',[q, q+1]))
                ZZ_instr_set.append(Gate('RZ', [q+1], angles=[psiZ]))
                ZZ_instr_set.append(Gate('CNOT',[q, q+1]))

            if self.h_ext != 0:
                P.add_instr(ext_instr_set)
            if self.JX !=0:
                P.add_instr(XX_instr_set)
            if self.JY !=0:
                P.add_instr(YY_instr_set)
            if self.JZ !=0:
                P.add_instr(ZZ_instr_set)
        return P

    def generate_local_circuits(self):

        ## Create circuits
        circuits = []
        for j in range(0, self.steps+1):
            print("Generating timestep {} circuit".format(j))
            with open(self.namevar,'a') as tempfile:
                tempfile.write("Generating timestep {} circuit\n".format(j))
            evolution_time = self.delta_t * j
            circuits.append(self.local_evolution_circuit(evolution_time))

        self.circuits_list=circuits


    def generate_ibm(self):
        #convert from local circuits to IBM-specific circuit
        #IBM imports 
        import qiskit as qk
        from qiskit.tools.monitor import job_monitor
        from qiskit.visualization import plot_histogram, plot_gate_map, plot_circuit_layout
        from qiskit import Aer, IBMQ, execute
        from qiskit.providers.aer import noise
        from qiskit.providers.aer.noise import NoiseModel
        from qiskit.circuit import quantumcircuit
        from qiskit.circuit import Instruction
        self.qr=qk.QuantumRegister(self.num_qubits, 'q')
        self.cr=qk.ClassicalRegister(self.num_qubits, 'c')
        if "y" in self.compile:
            ## Show available backends
            provider = qk.IBMQ.get_provider(group='open')
            provider.backends()

            #choose the device you would like to run on
            device = provider.get_backend(self.device_choice)

            #gather fidelity statistics on this device if you want to create a noise model for the simulator
            properties = device.properties()
            coupling_map = device.configuration().coupling_map

            #TO RUN ON THE SIMULATOR 
            #create a noise model to use for the qubits of the simulator
            noise_model = NoiseModel.from_backend(device)
            # Get the basis gates for the noise model
            basis_gates = noise_model.basis_gates

            # Select the QasmSimulator from the Aer provider
            simulator = Aer.get_backend('qasm_simulator')


            #To run on the quantum computer, assign a quantum computer of your choice as the backend 
            backend = provider.get_backend(self.device_choice)


        print("Creating IBM quantum circuit objects...")
        with open(self.namevar,'a') as tempfile:
            tempfile.write("Creating IBM quantum circuit objects...\n")
        name=0
        for circuit in self.circuits_list:
            propcirc = qk.QuantumCircuit(self.qr, self.cr)
            index=0
            for flip in self.flip_vec:
                if int(flip)==1:
                    propcirc.x(self.qr[index])
                    index+=1
                else: index+=1
            propcirc.barrier()
            for gate in circuit.gates:
                if "H" in gate.name:
                    propcirc.h(gate.qubits[0])
                elif "RZ" in gate.name:
                    propcirc.rz(gate.angles[0],gate.qubits[0])
                elif "RX" in gate.name:
                    propcirc.rx(gate.angles[0],gate.qubits[0])
                elif "CNOT" in gate.name:
                    propcirc.cx(gate.qubits[0],gate.qubits[1])
            propcirc.measure(self.qr,self.cr)
            self.ibm_circuits_list.append(propcirc)
        print("IBM quantum circuit objects created")
        with open(self.namevar,'a') as tempfile:
            tempfile.write("IBM quantum circuit objects created\n")

        if "y" in self.compile:
            if self.JZ != 0 and self.JX==self.JY==0 and self.h_ext!=0 and self.ext_dir=="X" and self.auto_smart_compile=="y":
                #TFIM
                print("TFIM detected, enabling DS compiler")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("TFIM detected, enabling DS compiler\n")
                temp=[]
                for circuit in self.ibm_circuits_list:
                    compiled=smart_compile(circuit,self.backend)
                    temp.append(compiled)
                self.ibm_circuits_list=temp

            elif self.default_compiler in "ds":
                temp=[]
                print("Compiling circuits...")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Compiling circuits...\n")
                for circuit in self.ibm_circuits_list:
                    compiled=smart_compile(circuit,self.backend)
                    temp.append(compiled)
                self.ibm_circuits_list=temp
                print("Circuits compiled successfully")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Circuits compiled successfully\n")
            elif self.default_compiler in "native":
                print("Transpiling circuits...")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Transpiling circuits...\n")
                temp=qk.compiler.transpile(self.ibm_circuits_list,backend=device,optimization_level=3)
                self.ibm_circuits_list=temp
                print("Circuits transpiled successfully")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Circuits transpiled successfully\n")


    def generate_rigetti(self):
        #Rigettti imports
        import pyquil
        from pyquil.quil import Program
        from pyquil.gates import H, RX, RZ, CZ, RESET, MEASURE
        from pyquil.api import get_qc

        print("Creating Pyquil program list...")
        with open(self.namevar,'a') as tempfile:
            tempfile.write("Creating Pyquil program list...\n")
        for circuit in self.circuits_list:
            p = pyquil.Program(RESET()) #compressed program
            ro = p.declare('ro', memory_type='BIT', memory_size=self.num_qubits)
            for gate in circuit.gates:
                if gate.name in "H":
                    p.inst(pyquil.gates.H(gate.qubits[0]))
                elif gate.name in "RZ":
                    p.inst(pyquil.gates.RZ(gate.angles[0],gate.qubits[0]))
                elif gate.name in "RX":
                    p.inst(pyquil.gates.RX(gate.angles[0],gate.qubits[0]))
                elif gate.name in "CNOT":
                    p.inst(pyquil.gates.CNOT(gate.qubits[0],gate.qubits[1]))
            for i in range(self.num_qubits):
                p.inst(pyquil.gates.MEASURE(i,ro[i]))
            p.wrap_in_numshots_loop(self.shots)
            self.rigetti_circuits_list.append(p)

        if "y" in self.compile:
            qc=get_qc(self.device_choice)
            if self.JZ != 0 and self.JX==self.JY==0 and self.h_ext!=0 and self.ext_dir=="X" and self.auto_smart_compile=="y":
                #TFIM
                print("TFIM detected, enabling DS compiler")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("TFIM detected, enabling DS compiler\n")
                temp=[]
                for circuit in self.rigetti_circuits_list:
                    temp.append(smart_compile(circuit,self.backend,self.shots))
                self.rigetti_circuits_list=temp

            elif self.default_compiler in "ds":
                temp=[]
                print("Compiling circuits...")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Compiling circuits...\n")
                for circuit in self.rigetti_circuits_list:
                    temp.append(smart_compile(circuit,self.backend,self.shots))
                self.rigetti_circuits_list=temp
                print("Circuits compiled successfully")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Circuits compiled successfully\n")
            elif self.default_compiler in "native":
                temp=[]
                print("Transpiling circuits...")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Transpiling circuits...\n")
                for circuit in self.rigetti_circuits_list:
                    circ = qc.compile(circuit)
                    temp.append(circ)
                self.rigetti_circuits_list=temp
                print("Circuits transpiled successfully")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Circuits transpiled successfully\n")

        print("Pyquil program list created successfully")
        with open(self.namevar,'a') as tempfile:
            tempfile.write("Pyquil program list created successfully\n")

    def generate_cirq(self):
        #Cirq imports
        import cirq
        print("Creating Cirq circuit list...")
        with open(self.namevar,'a') as tempfile:
            tempfile.write("Creating Cirq circuit list...\n")
        for circuit in self.circuits_list:
            c=cirq.Circuit()
            qubit_list=cirq.LineQubit.range(self.num_qubits)
            gate_list=[]
            for gate in circuit.gates:
                if gate.name in "H":
                    gate_list.append(cirq.H(qubit_list[gate.qubits[0]]))
                elif gate.name in "RZ":
                    gate_list.append(cirq.rz(gate.angles[0])(qubit_list[gate.qubits[0]]))
                elif gate.name in "RX":
                    gate_list.append(cirq.rx(gate.angles[0])(qubit_list[gate.qubits[0]]))
                elif gate.name in "CNOT":
                    gate_list.append(cirq.CNOT(qubit_list[gate.qubits[0]],qubit_list[gate.qubits[1]]))
            for i in range(self.num_qubits):
                gate_list.append(cirq.measure(qubit_list[i]))
            c.append(gate_list,strategy=cirq.InsertStrategy.EARLIEST)
            self.cirq_circuits_list.append(c)
        print("Successfully created Cirq circuit list")
        with open(self.namevar,'a') as tempfile:
            tempfile.write("Successfully created Cirq circuit list\n")

    def generate_circuits(self):
        # self.imports()
        if len(self.circuits_list)==0:
            self.generate_local_circuits()

        if self.backend in "ibm":
            self.generate_ibm()
        if self.backend in "rigetti":
            self.generate_rigetti()
        if self.backend in "cirq":
            self.generate_cirq()

    def connect_IBM(self,api_key=None, overwrite=False):
        import qiskit as qk
        if api_key != None:
            if overwrite==False:
                qk.IBMQ.save_account(api_key) ## only run once!
            else:
                qk.IBMQ.save_account(api_key,overwrite=True) ## only run once!
     #qk.IBMQ.delete_accounts() ## only run if you need to use a new token
        qk.IBMQ.load_account()


    def parameters(self):
        print("Current model parameters:\n\nH_BAR = {}\nJX = {}\nJY = {}\nJZ = {}\nh_ext = {}\next_dir = {}".format(self.H_BAR,self.JX,self.JY,self.JZ,self.h_ext,self.ext_dir))
        print("num_qubits = {}\ninitial_spins = {}\ndelta_t = {}\nsteps = {}\nQCQS = {}\nshots = {}\nnoise_choice = {}".format(self.num_qubits,self.initial_spins,self.delta_t,self.steps,self.QCQS,self.shots,self.noise_choice))
        print("device choice = {}\nplot_flag = {}\nfreq = {}\ntime_dep_flag = {}\ncustom_time_dep = {}\n".format(self.device_choice,self.plot_flag,self.freq,self.time_dep_flag,self.custom_time_dep)) 
        print("compile = {}\nauto_smart_compile = {}\ndefault_compiler = {}".format(self.compile,self.auto_smart_compile,self.default_compiler))
        #this is missing some of the latest parameter additions

    def results(self):
        return self.result_matrix

    def return_circuits(self):
        if self.backend in "ibm":
            if len(self.ibm_circuits_list)==0:
                self.generate_circuits()
            return self.ibm_circuits_list
        elif self.backend in "rigetti":
            if len(self.rigetti_circuits_list)==0:
                self.generate_circuits()
            return self.rigetti_circuits_list
        elif self.backend in "cirq":
            if len(self.cirq_circuits_list)==0:
                self.generate_circuits()
            return self.cirq_circuits_list

    def average_magnetization(self,result: dict, shots: int, qub: int):
      """Compute average magnetization from results of qk.execution.
      Args:
      - result (dict): a dictionary with the counts for each qubit, see qk.result.result module
      - shots (int): number of trials
      Return:
      - average_mag (float)
      """
      mag = 0
      for spin_str, count in result.items():
        spin_int = [1 - 2 * float(spin_str[qub])]
        #print(spin_str)
        mag += (sum(spin_int) / len(spin_int)) * count
      average_mag = mag / shots
      return average_mag














    def run_circuits(self):
        if "y" in self.plot_flag:
            import matplotlib.pyplot as plt
        if self.backend in "ibm":
            import qiskit as qk
            from qiskit.tools.monitor import job_monitor
            from qiskit.visualization import plot_histogram, plot_gate_map, plot_circuit_layout
            from qiskit import Aer, IBMQ, execute
            from qiskit.providers.aer import noise
            from qiskit.providers.aer.noise import NoiseModel
            from qiskit.circuit import quantumcircuit
            from qiskit.circuit import Instruction
            ## Show available backends
            provider = qk.IBMQ.get_provider(group='open')
            provider.backends()

            #choose the device you would like to run on
            device = provider.get_backend(self.device_choice)

            #gather fidelity statistics on this device if you want to create a noise model for the simulator
            properties = device.properties()
            coupling_map = device.configuration().coupling_map

            #TO RUN ON THE SIMULATOR 
            #create a noise model to use for the qubits of the simulator
            noise_model = NoiseModel.from_backend(device)
            # Get the basis gates for the noise model
            basis_gates = noise_model.basis_gates

            # Select the QasmSimulator from the Aer provider
            simulator = Aer.get_backend('qasm_simulator')


            #To run on the quantum computer, assign a quantum computer of your choice as the backend 
            backend = provider.get_backend(self.device_choice)

            #CHOOSE TO RUN ON QUANTUM COMPUTER OR SIMULATOR
            if self.QCQS in ["QC"]:
                #quantum computer execution
                job = qk.execute(self.ibm_circuits_list, backend=backend, shots=self.shots)
                job_monitor(job)

            elif self.QCQS in ["QS"]:
                #simulator execution
                if self.noise_choice in ["y"]:
                    print("Running noisy simulator job...")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Running noisy simulator job...\n")
                    result_noise = execute(self.ibm_circuits_list, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates,shots=self.shots).result()
                    print("Noisy simulator job successful")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Noisy simulator job successful\n")
                elif self.noise_choice in ["n"]:
                    print("Running noiseless simulator job...")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Running noiseless simulator job...\n")
                    result_noise=execute(self.ibm_circuits_list,simulator,coupling_map=coupling_map,basis_gates=basis_gates,shots=self.shots).result()
                    print("Noiseless simulator job successful")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Noiseless simulator job successful\n")
                else: 
                    print("Please enter either y or n for the simulator noise query")
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Please enter either y or n for the simulator noise query\n")
            else:
                print("Please enter either QC or QS")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Please enter either QC or QS\n")


            
                


            #Post Processing Depending on Choice
            self.result_out_list=[]
            if self.QCQS in ["QS"]:
                #SIMULATOR POST PROCESSING
                for j in range(self.num_qubits):
                    avg_mag_sim = []
                    temp = []
                    i = 1
                    print("Post-processing qubit {} data".format(j+1))
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Post-processing qubit {} data\n".format(j+1))
                    for c in self.ibm_circuits_list:
                        result_dict = result_noise.get_counts(c)
                        temp.append(self.average_magnetization(result_dict, self.shots,j))
                        if i % (self.steps+1) == 0:
                            avg_mag_sim.append(temp)
                            temp = []
                        i += 1
                    # time_vec=np.linspace(0,total_t,steps)
                    # time_vec=time_vec*JX/H_BAR
                    if "y" in self.plot_flag:
                        plt.figure()
                        plt.plot(range(self.steps+1), avg_mag_sim[0])
                        plt.xlabel("Simulation Timestep")
                        plt.ylabel("Average Magnetization")
                        plt.savefig("Data/Simulator_result_qubit{}.png".format(j+1))
                        plt.close()
                    self.result_out_list.append(avg_mag_sim[0])
                    np.savetxt("Data/Qubit {} Average Magnetization Data.txt".format(j+1),avg_mag_sim[0])
                self.result_matrix=np.stack(self.result_out_list)
                print("Done")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Done\n")

            elif self.QCQS in ["QC"]:
                #QUANTUM COMPUTER POST PROCESSING
                for j in range(self.num_qubits):
                    results = job.result()        
                    avg_mag_qc = []
                    temp = []
                    i = 1
                    print("Post-processing qubit {} data".format(j+1))
                    with open(self.namevar,'a') as tempfile:
                        tempfile.write("Post-processing qubit {} data\n".format(j+1))
                    for c in self.ibm_circuits_list:
                            result_dict = results.get_counts(c)
                            temp.append(self.average_magnetization(result_dict, self.shots,j))
                            if i % (self.steps+1) == 0:
                                    avg_mag_qc.append(temp)
                                    temp = []
                            i += 1
                    
                    # QC
                    if "y" in self.plot_flag:
                        plt.figure()
                        plt.plot(range(self.steps+1), avg_mag_qc[0])
                        plt.xlabel("Simulation Timestep")
                        plt.ylabel("Average Magnetization")
                        plt.savefig("Data/QC_result_qubit{}.png".format(j+1))
                        plt.close()
                    self.result_out_list.append(avg_mag_qc[0])
                    np.savetxt("Data/Qubit {} Average Magnetization Data.txt".format(j+1),avg_mag_qc[0])
                self.result_matrix=np.stack(self.result_out_list)           
                print("Done")
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Done\n")
        elif "rigetti" in self.backend:
            print("Running Pyquil programs...")
            with open(self.namevar,'a') as tempfile:
                tempfile.write("Running Pyquil programs...\n")
            qc=get_qc(self.device_choice)
            results_list=[]
            first_ind=0
            #each circuit represents one timestep
            for circuit in self.rigetti_circuits_list:
                # print("Ay I got a circuit here")
                temp=qc.run(circuit)
                results_list.append(temp)

            for i in range(self.num_qubits):
                print("Post-processing qubit {} data...".format(i+1))
                with open(self.namevar,'a') as tempfile:
                    tempfile.write("Post-processing qubit {} data...\n".format(i+1))
                qubit_specific_row=np.zeros(len(results_list))
                for j in range(len(self.rigetti_circuits_list)):
                    results=results_list[j]

                    summation=0
                    for array in results:
                        summation+=(1-2*array[i])

                    summation=summation/len(results) #average over the number of shots

                    qubit_specific_row[j]=summation
                if first_ind==0:
                    self.result_matrix=qubit_specific_row
                    first_ind+=1
                else:
                    self.result_matrix=np.vstack((self.result_matrix,qubit_specific_row))
                if "y" in self.plot_flag:
                    plt.figure()
                    xaxis=np.linspace(0,self.steps,num=self.steps+1)
                    plt.plot(qubit_specific_row)
                    plt.xlabel("Simulation Timestep")
                    plt.ylabel("Average Magnetization")
                    plt.savefig("Data/Result_qubit{}.png".format(i+1))
                    plt.close()
            print("Done")
            with open(self.namevar,'a') as tempfile:
                tempfile.write("Done\n")
