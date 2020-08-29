from quantum_circuits import Program
import numpy as np
def ds_compile(circ_obj, circ_type, shots=1):
    if (circ_type == "ibm"):
        return ds_compile_ibm(circ_obj,shots=1)
    if (circ_type == "rigetti"):
        return ds_compile_rigetti(circ_obj,shots=1)
    else:
        print("invalid circuit type. Use rigetti or ibm")

def ds_compile_ibm(circ_obj,shots=1):
    #IBM imports 
    import qiskit as qk
    from qiskit.circuit import quantumcircuit
    nqubits=circ_obj.num_qubits    
    #Read the gate in right vector form
    # G = Gate type
    # TH = Angle of rotation ! if no angle rotation then TH = 0
    # TH2 = 2nd angle of rotation (used in U2 and U3 gates)
    # TH3 = 3rd angle of rotation (used in U3 gates)
    # AC1 = qubit on which action is happening
    # AC2 = qubit on which controlled action is happening
    instr_list=circ_obj.data
    count = len(instr_list)

    G = ["" for x in range(count)]
    G = list(G)
    AC1 = np.zeros(shape=(count),dtype=np.int) 
    AC2 = np.zeros(shape=(count),dtype=np.int) 
    TH = np.zeros(shape=(count))
    i = 0
    for instr in instr_list:
        G[i] = 0
        TH[i] = 0
        AC1[i] = 0
        AC2[i] = 0
        name = instr[0].name
        if name == "h":
            G[i]="H"
            TH[i] = 0
            AC1[i] = instr[1][0].index
            AC2[i] = 0
        if name == "rz":
            G[i] = "RZ"
            TH[i] = instr[0].params[0]
            AC1[i] = instr[1][0].index
            AC2[i] = 0
        if name == "cx":
            G[i] = "CNOT"
            TH[i] = 0
            AC1[i] = instr[1][0].index
            AC2[i] = instr[1][1].index
        if name == "measure":
            G[i] = "MEASURE"
            TH[i] = 0
            AC1[i] = 0
            AC2[i] =0
        if name == "rx":
            G[i] = "RX"
            TH[i] = instr[0].params[0]
            AC1[i] = instr[1][0].index
            AC2[i] = 0
        i = i+1


    #Use RX = H RZ H
    i=0
    while G[i]!= "MEASURE":
        if G[i]=="RX":
            G[i]="H"
            intermed_angle=float(TH[i].copy())
            intermed_qubit=int(AC1[i].copy())
            G.insert(i,"RZ")
            TH=np.insert(TH,i,intermed_angle)
            AC1=np.insert(AC1,i,intermed_qubit)
            AC2=np.insert(AC2,i,0)
            G.insert(i,"H")
            TH=np.insert(TH,i,0)
            AC1=np.insert(AC1,i,intermed_qubit)
            AC2=np.insert(AC2,i,0)
        i=i+1
            
    #Omit last and second-to-last CNOT for each qubit
    for qub in range(0,nqubits+1):
        i=-1
        count=0
        while count<=1 and i>=-int(len(G)):
            if G[i] == "CNOT" and AC1[i]==qub and AC2[i]==qub+1:
                del G[i]
                TH=np.delete(TH,i)
                AC1=np.delete(AC1,i)
                AC2=np.delete(AC2,i)
                count=count+1
            i=i-1

    #Omit last RZ for each qubit
    for qub in range(0,nqubits+1):
        i=-1
        while i>=-int(len(G)):
            if G[i] == "H" and AC1[i]==qub:
                break
            if G[i] == "RZ" and AC1[i]==qub:
                G[i]  = "NULL"
                break
            i=i-1            



    #Use CNOT (0,1) ->  H(0) H(1) CNOT(1,0) H(0) H(1)
    i=0
    while G[i] != "MEASURE":
        if G[i]=="CNOT" and (G[i+1]=="H" and G[i+2]=="H" and AC1[i+1]==AC1[i] and AC1[i+2]==AC2[i])==False:
            G[i]="H"
            flag1=int(AC1[i])
            flag2=int(AC2[i])
            AC2[i]=0
            G.insert(i,"H")
            TH=np.insert(TH,i,0)
            AC1=np.insert(AC1,i,flag2)
            AC2=np.insert(AC2,i,0)
            G.insert(i,"CNOT")
            TH=np.insert(TH,i,0)
            AC1=np.insert(AC1,i,flag2)
            AC2=np.insert(AC2,i,flag1)
            G.insert(i,"H")
            TH=np.insert(TH,i,0)
            AC1=np.insert(AC1,i,flag1)
            AC2=np.insert(AC2,i,0)
            G.insert(i,"H")
            TH=np.insert(TH,i,0)
            AC1=np.insert(AC1,i,flag2)
            AC2=np.insert(AC2,i,0)
        i=i+1

    #Rearrange circuits to put successive Hadamard gates in order
    i=0
    while G[i] != "MEASURE":
        if G[i]=="H":
            flag=AC1[i]
            j=i+1
            boolean=0
            while G[j] != "MEASURE" and boolean ==0:
                if AC1[j]==flag and G[j] == "H":
                    boolean=1
                    del G[j]
                    TH=np.delete(TH,j)
                    AC1=np.delete(AC1,j)
                    AC2=np.delete(AC2,j)
                    G.insert(i,"H")
                    TH=np.insert(TH,i,0)
                    AC1=np.insert(AC1,i,flag)
                    AC2=np.insert(AC2,i,0)
                if AC1[j]==flag and G[j] != "H":
                    break
                j=j+1
        i=i+1


    #Use successive Hadamard annihilation
    i=0
    while G[i]!= "MEASURE":
        if G[i]=="H" and G[i+1] == "H" and AC1[i]==AC1[i+1]:
            del G[i]
            TH=np.delete(TH,i)
            AC1=np.delete(AC1,i)
            AC2=np.delete(AC2,i)
            del G[i]
            TH=np.delete(TH,i)
            AC1=np.delete(AC1,i)
            AC2=np.delete(AC2,i)
            i=i-1
        i=i+1


    #Convert HRZ(theta)H to RZ(pi/2)RX(pi/2)RZ(theta+pi)RX(pi/2)RZ(pi/2)
    i=0
    while G[i] != "MEASURE":
        if (G[i] == "H" and G[i+1] == "RZ" and G[i+2]=="H" and AC1[i] == AC1[i+1] and AC1[i+1]== AC1[i+2]):
            theta = TH[i+1]
            q = AC1[i]
            G[i]="RZ"
            TH[i]=1.57079632679
            del G[i+1]
            TH=np.delete(TH,i+1)
            AC1=np.delete(AC1,i+1)
            AC2=np.delete(AC2,i+1)
            del G[i+1]
            TH=np.delete(TH,i+1)
            AC1=np.delete(AC1,i+1)
            AC2=np.delete(AC2,i+1)
            G.insert(i,"RX")
            TH=np.insert(TH,i,1.57079632679)
            AC1=np.insert(AC1,i,q)
            AC2=np.insert(AC2,i,q)
            G.insert(i,"RZ")
            TH=np.insert(TH,i,theta+(2.0*1.57079632679))
            AC1=np.insert(AC1,i,q)
            AC2=np.insert(AC2,i,q)
            G.insert(i,"RX")
            TH=np.insert(TH,i,1.57079632679)
            AC1=np.insert(AC1,i,q)
            AC2=np.insert(AC2,i,q)
            G.insert(i,"RZ")
            TH=np.insert(TH,i,1.57079632679)
            AC1=np.insert(AC1,i,q)
            AC2=np.insert(AC2,i,q)

            #move leftmost RZ of set across control bit if possible
            for j in range(i-1,0,-1):
                if AC1[j] == AC1[i]:
                    if G[j] == "CNOT":
                        for k in range(j-1,0,-1):
                            if AC1[k] == AC1[i]:
                                if G[k] == "RZ":
                                    TH[k]=TH[k]+TH[i]
                                    del G[i]
                                    TH=np.delete(TH,i)
                                    AC1=np.delete(AC1,i)
                                    AC2=np.delete(AC2,i)
                                else: break
                    else: break

            #move rightmost RZ of set across control bit if possible
            for j in range(i+4,len(G)):
                if AC1[j] == AC1[i+3]:
                    if G[j] == "CNOT":
                        for k in range(j+1,len(G)):
                            if AC1[k] == AC1[i+3]:
                                if G[k] == "RZ":
                                    TH[k]=TH[k]+TH[i+3]
                                    del G[i+3]
                                    TH=np.delete(TH,i+3)
                                    AC1=np.delete(AC1,i+3)
                                    AC2=np.delete(AC2,i+3)
                                else: break
                            if AC2[k] == AC1[i+3]:
                                break
                    else: break


        i=i+1

    #convert remaining HRZ or H to native gates
    i=0
    while G[i] != "MEASURE":
        if G[i]=="H":
            q = AC1[i]
            j=i+1
            flag = 1
            while G[j] != "MEASURE":
                if AC1[j] == AC1[i]:
                    #change HRZ to native gates
                    if G[j]=="RZ":
                        G[i] = "RZ"
                        theta = TH[j]
                        TH[i]=1.57079632679
                        del G[j]
                        TH=np.delete(TH,j)
                        AC1=np.delete(AC1,j)
                        AC2=np.delete(AC2,j)
                        G.insert(i+1,"RX")
                        TH=np.insert(TH,i+1,1.57079632679)
                        AC1=np.insert(AC1,i+1,q)
                        AC2=np.insert(AC2,i+1,0)
                        G.insert(i+2,"RZ")
                        TH=np.insert(TH,i+2,theta+1.57079632679)
                        AC1=np.insert(AC1,i+2,q)
                        AC2=np.insert(AC2,i+2,0)
                        flag = 0
                        break
                    else: break
                j=j+1
            #change H to native gates    
            if (flag):
                G[i] = "RZ"
                TH[i]=1.57079632679
                G.insert(i+1,"RX")
                TH=np.insert(TH,i+1,1.57079632679)
                AC1=np.insert(AC1,i+1,q)
                AC2=np.insert(AC2,i+1,0)
                G.insert(i+2,"RZ")
                TH=np.insert(TH,i+2,1.57079632679)
                AC1=np.insert(AC1,i+2,q)
                AC2=np.insert(AC2,i+2,0)  
            #compress successive RZs
            if (G[i-1] == "RZ" and AC1[i-1] == AC1[i]):
                TH[i-1] = TH[i-1]+TH[i]
                del G[i]
                TH=np.delete(TH,i)
                AC1=np.delete(AC1,i)
                AC2=np.delete(AC2,i)
            #if (G[i+3] == "RZ"):
            #    TH[i+2] = TH[i+2]+TH[i+3]
            #    del G[i+3]
            #    TH=np.delete(TH,i+3)
            #    AC1=np.delete(AC1,i+3)
            #    AC2=np.delete(AC2,i+3)

        i=i+1

    #Omit first RZs
    for qub in range(0,nqubits):
        i=0
        while G[i] != "MEASURE":
            if G[i]=="RZ" and AC1[i]==qub:
                del G[i]
                TH=np.delete(TH,i)
                AC1=np.delete(AC1,i)
                AC2=np.delete(AC2,i)
            if (G[i]=="RX" and AC1[i]==qub) or (G[i]=="CNOT" and (AC1[i]==qub or AC2[i]==qub)):
                break
            i=i+1

    #Omit last RZ for each qubit
    for qub in range(0,nqubits+1):
        i=-1
        while i>=-int(len(G)):
            if G[i] == "H" and AC1[i]==qub:
                break
            if G[i] == "RZ" and AC1[i]==qub:
                G[i]  = "NULL"
                break
            i=i-1  

    #build output circuit
    qr = qk.QuantumRegister(nqubits, 'q')
    cr = qk.ClassicalRegister(nqubits, 'c')

    circuit = qk.QuantumCircuit(qr, cr)

    for i in range(len(G)):
        if (G[i] == "RX"):
            circuit.rx(TH[i], int(AC1[i]))
        if (G[i] == "RZ"):
            circuit.rz(TH[i], int(AC1[i]))
        if (G[i] == "CNOT"):
            circuit.cx(int(AC1[i]), int(AC2[i]))
        if (G[i] == "H"):
            circuit.h(int(AC1[i]))


    circuit.measure(qr, cr)

    return circuit


def ds_compile_rigetti(circ_obj,shots=1):
    import pyquil
    from pyquil.gates import RX, RZ, CZ, MEASURE
    from pyquil.api import get_qc
    nqubits=len(circ_obj.get_qubits())
    lineList = [str(instr) for instr in circ_obj]
    count = len(lineList)
    #Read the gate in right vector form
    # G = Gate type
    # TH = Angle of rotation ! if no angle rotation then TH = 0
    # AC1 = qubit on which action is happening
    # AC2 = qubit on which controlled action is happening

    G = ["" for x in range(count)]
    G = list(G)
    AC1 = np.zeros(shape=(count),dtype=np.int) 
    AC2 = np.zeros(shape=(count),dtype=np.int) 
    TH = np.zeros(shape=(count))
    for i in range (0,count):
      G[i] = 0
      TH[i] = 0
      AC1[i] = 0
      AC2[i] = 0
      if lineList[i][0:1] == "H":
        G[i]="H"
        TH[i] = 0
        AC1[i] = lineList[i][2:3]
        AC2[i] = 0
      if lineList[i][0:2] == "RZ":
        G[i] = "RZ"
        TH[i] = lineList[i][lineList[i].find("(")+1:lineList[i].find(")")]
        AC1[i] = lineList[i][-1]
        AC2[i] = 0
      if lineList[i][0:2] == "RX":
        G[i] = "RX"
        TH[i] = lineList[i][lineList[i].find("(")+1:lineList[i].find(")")]
        AC1[i] = lineList[i][-1]
        AC2[i] = 0
      if lineList[i][0:4] == "CNOT":
        G[i] = "CNOT"
        TH[i] = 0
        AC1[i] = lineList[i][5:6] 
        AC2[i] = lineList[i][7:8]
      if lineList[i][0:7] == "MEASURE":
        G[i] = "MEASURE"
        TH[i] = 0
        AC1[i] = 0
        AC2[i] =0

    #Use RX = H RZ H
    i=0
    while G[i]!= "MEASURE":
        if G[i]=="RX":
            G[i]="H"
            intermed_angle=TH[i].copy()
            intermed_qubit=AC1[i].copy()
            G.insert(i,"RZ")
            TH=np.insert(TH,i,intermed_angle)
            AC1=np.insert(AC1,i,intermed_qubit)
            AC2=np.insert(AC2,i,0)
            G.insert(i,"H")
            TH=np.insert(TH,i,0)
            AC1=np.insert(AC1,i,intermed_qubit)
            AC2=np.insert(AC2,i,0)
        i=i+1


    # Use CNOT = H CZ H 
    i = 0
    while G[i] != "MEASURE":
      if G[i] == "CNOT":
         G[i] = "CZ"
         G.insert(i+1,"H")
         TH = np.insert(TH,i+1,0) 
         AC1 = np.insert(AC1,i+1,AC2[i]) 
         AC2 = np.insert(AC2,i+1,0)
         G.insert(i,"H")
         TH = np.insert(TH,i,0) 
         AC1 = np.insert(AC1,i,AC2[i]) 
         AC2 = np.insert(AC2,i,0)
      i = i+1  

    # Last and second last CNOT can be ommited  
    maxq = max(max(AC1),max(AC2))
    remember = np.zeros(shape=(2,maxq),dtype=np.int)
    for mm in range (0,maxq+1):
     i = 0
     while G[i] != "MEASURE":
          if G[i] == "CZ" and AC1[i] == mm and AC2[i] == mm+1:
                j = i+1
                while G[j] != "MEASURE":
                      if G[j] == "CZ" and AC1[j] == mm and AC2[j] == mm+1:
                           remember[0][mm] = i; remember[1][mm] = j;
                      j = j+1
          i = i+1 

    for nn in range (maxq-1,-1,-1):
     for mm in range (1,-1,-1):
    #   print(mm,nn)
      del G[remember[mm][nn]];TH = np.delete(TH,remember[mm][nn]);
      AC1 = np.delete(AC1,remember[mm][nn]); AC2 = np.delete(AC2,remember[mm][nn])


    # Use H*H = I but make sure it can only happen if no gate is 
    # present in between 
    i = 0
    while G[i] != "MEASURE":
      if G[i] == "H":
        flag = 0
        #print(G[i],TH[i],AC1[i],AC2[i],"before start")
        j = i+1
        while G[j] != "MEASURE":
           if ((G[j] == "CZ" and AC1[j] == AC1[i]) or (G[j] == "CZ" and AC2[j] == AC1[i]) or (G[j] == "RZ" and AC1[j] == AC1[i])) :
              break
           if G[j] == G[i] and AC1[j] == AC1[i] :
              #print(G[i],TH[i],AC1[i],AC2[i],"before")
              del G[j]
              TH = np.delete(TH,j)
              AC1 = np.delete(AC1,j)
              AC2 = np.delete(AC2,j)
              #print(G[i],TH[i],AC1[i],AC2[i],"after")
              del G[i]
              TH = np.delete(TH,i)
              AC1 = np.delete(AC1,i)
              AC2 = np.delete(AC2,i)
              flag = 2
           j = j+1
           if flag ==2:
              break 
      i = i + 1


    # Use CZ H RZ H CZ = RZ(pi/2) CZ RX(pi/2) RZ RX(-pi2) CZ RZ(-pi/2)
    i = 0 
    while G[i] != "MEASURE":
        if (G[i] == "CZ" and G[i+1] == "H" and AC2[i] == AC1[i+1] and G[i+2] == "RZ" and AC2[i] == AC1[i+2] and G[i+3] == "H" and AC2[i] == AC1[i+3] and G[i+4] == "CZ" and AC2[i] == AC2[i+4]):
              G[i+1] = "RX"; TH[i+1] = 1.57079632679; 
              G[i+3] = "RX"; TH[i+3] = -1.57079632679;
              G.insert(i+5,"RZ"); TH = np.insert(TH,i+5,-1.57079632679); 
              AC1 = np.insert(AC1,i+5,AC2[i]); AC2 = np.insert(AC2,i+5,0);
              G.insert(i,"RZ"); TH = np.insert(TH,i,1.57079632679); 
              AC1 = np.insert(AC1,i,AC2[i]); AC2 = np.insert(AC2,i,0);
              # print("loop activated")
        i = i+1


    # Use H = RZ(pi/2) RX(pi/2) RZ(pi/2)
    i = 0
    while G[i] !="MEASURE":
        if (G[i] == "H"):
              flag = AC1[i]
              G[i] = "RZ"; TH[i] = 1.57079632679 ;
              G.insert(i,"RX");TH = np.insert(TH,i,1.57079632679);
              AC1 = np.insert(AC1,i,flag); AC2 = np.insert(AC2,i,0); 
              G.insert(i,"RZ");TH = np.insert(TH,i,1.57079632679); 
              AC1 = np.insert(AC1,i,flag); AC2 = np.insert(AC2,i,0); 
        i = i+1 


    # Compress RZ gates
    loop_flag = 0
    for mm in range (0,1000):
     i = 0 
     while G[i] !="MEASURE": 
         if (G[i] == "RZ"):
               j = i+1
               flag = 0
               #print(flag,"flag")
               while G[j] !="MEASURE":
                     if (G[j] == "RX" and AC1[j] == AC1[i]):
                          flag = 2
                     if (G[j] == "RZ" and AC1[j] == AC1[i]):
                          TH[i] = TH[i]+TH[j]; 
                          del G[j];TH = np.delete(TH,j);
                          AC1 = np.delete(AC1,j); AC2 = np.delete(AC2,j) 
                          flag = 2
                          loop_flag = 3
                     j = j+1
                     if(flag == 2):
                          break
         if (G[i] == "RZ" and TH[i]== 0.0):
               del G[i];TH = np.delete(TH,i);
               AC1 = np.delete(AC1,i); AC2 = np.delete(AC2,i)
         i = i +1
     if(loop_flag == 0):
         break  
     if(mm ==1000 and loop_flag==3):
         print("more RZ compression are left be carefull!!")


    i = 0
    while G[i] != "MEASURE":
        if (G[i] == "RX" and TH[i] == 1.57079632679):
            i1 = i+1
            while G[i1] != "MEASURE":
                if (G[i1] == "RX" and AC1[i1] == AC1[i] or G[i1] == "CZ" and AC1[i1] == AC1[i] or G[i1] == "RZ" and TH[i1] != 3.14159265358 and AC1[i1] == AC1[i] or G[i1] == "CZ" and AC2[i1] == AC1[i]):
                    break
                if (G[i1] == "RZ" and TH[i1] == 3.14159265358 and AC1[i1] == AC1[i]):
                    i2 = i1+1
                    while G[i2] != "MEASURE":
                        if (G[i2] == "RX" and AC1[i2] == AC1[i] or G[i2] == "RZ" and AC1[i2] == AC1[i] ):
                            break
                        if (G[i2] == "CZ" and AC1[i2] == AC1[i] and G[i2+4] == "CZ" and AC1[i2+4] == AC1[i]):
                            i3 = i2 +5
                            while G[i3] != "MEASURE":
                                if(G[i3] == "RZ" and AC1[i3] == AC1[i]+1 or G[i3] == "CZ" and AC1[i3] == AC1[i]+1 or G[i3] == "RX" and TH[i3] != 1.57079632679 and AC1[i3] == AC1[i]+1 ):
                                    break
                                if(G[i3] == "RX" and TH[i3] == 1.57079632679 and AC1[i3] == AC1[i]+1) :
                                    i4 = i2 + 5
                                    while G[i4] != "MEASURE":
                                        if (G[i4] == "RZ" and AC1[i4] == AC1[i] or G[i4] == "CZ" and AC1[i4] == AC1[i] or G[i4] == "RX" and TH[i4] != 1.57079632679 and AC1[i4] == AC1[i] ):
                                            break
                                        if(G[i4] == "RX" and TH[i4] == 1.57079632679 and AC1[i4] == AC1[i]) :
                                            AC1[i2+1] = AC1[i];AC1[i2+2] = AC1[i];AC1[i2+3] = AC1[i]
                                            G[i4] = "RZ"; TH[i4] = 3.14159265358;
                                            del G[i3];TH = np.delete(TH,i3)
                                            AC1 = np.delete(AC1,i3); AC2 = np.delete(AC2,i3)
                                            G.insert(i2,"RX")
                                            TH = np.insert(TH,i2,1.57079632679)           
                                            AC1 = np.insert(AC1,i2,AC2[i2]); AC2 = np.insert(AC2,i2,0)
                                            del G[i1];TH = np.delete(TH,i1)
                                            AC1 = np.delete(AC1,i1); AC2 = np.delete(AC2,i1)
                                            del G[i];TH = np.delete(TH,i)
                                            AC1 = np.delete(AC1,i); AC2 = np.delete(AC2,i)
                                            break
                                        i4 = i4 +1
                                i3 = i3 + 1
                        i2 = i2 + 1
                i1 = i1 + 1 
        i = i+1


    # Compress RZ gates                                                            
    loop_flag = 0                                                                  
    for mm in range (0,1000):                                                      
        i = 0                                                                         
        while G[i] !="MEASURE":                                                       
            if (G[i] == "RZ"):                                                        
                j = i+1                                                             
                flag = 0                                                            
                while G[j] !="MEASURE":                                             
                    if (G[j] == "RX" and AC1[j] == AC1[i]):                       
                        flag = 2                                                 
                    if (G[j] == "RZ" and AC1[j] == AC1[i]):                       
                        TH[i] = TH[i]+TH[j];                                     
                        del G[j];TH = np.delete(TH,j);                           
                        AC1 = np.delete(AC1,j); AC2 = np.delete(AC2,j)           
                        flag = 2                                                 
                        loop_flag = 3                                            
                    j = j+1                                                       
                    if(flag == 2):                                                
                        break                                                    
            if (G[i] == "RZ" and TH[i]== 0.0):                                        
                del G[i];TH = np.delete(TH,i);                                      
                AC1 = np.delete(AC1,i); AC2 = np.delete(AC2,i)                      
            i = i +1                                                                  
        if(loop_flag == 0):                                                           
            break                                                                     
        if(mm ==1000 and loop_flag==3):                                               
            print("more RZ compression are left be carefull!!")


    # Use RZ(theta) RX RZ(pi) = RZ(theta-pi) RX(-pi2)
    i = 0
    while G[i] != "MEASURE":
        if (G[i] == "RZ"):
            loop_breaker = 0
            i1 = i+1
            while G[i1] != "MEASURE":
                if (G[i1] == "RX" and TH[i1] != 1.57079632679 and AC1[i1] == AC1[i]):
                    break
                if (G[i1] == "RX" and TH[i1] == 1.57079632679 and AC1[i1] == AC1[i]):
                    i2 = i1 + 1
                    while G[i2] != "MEASURE":
                        if (G[i2] == "RZ" and TH[i2] == 3.14159265358 and AC1[i2] == AC1[i]):
                            TH[i] = TH[i]+TH[i2]
                            TH[i1] = -TH[i1]
                            del G[i2];TH = np.delete(TH,i2);
                            AC1 = np.delete(AC1,i2); AC2 = np.delete(AC2,i2);
                            loop_breaker = 3
                            break
                        elif (G[i2] == "RZ" and TH[i2] != 3.14159265358 and AC1[i2] == AC1[i]):
                            loop_breaker = 3
                            break
                        i2 = i2 + 1                         
                if (loop_breaker == 3):
                    break
                i1 = i1 + 1
        i = i + 1


    p = pyquil.Program() #compressed program
    ro = p.declare('ro', memory_type='BIT', memory_size=nqubits)

    for i in range(len(G)):
        if (G[i] == "RX"):
            p.inst(pyquil.gates.RX(TH[i], int(AC1[i])))
        if (G[i] == "RZ"):
            p.inst(pyquil.gates.RZ(TH[i], int(AC1[i])))
        if (G[i] == "CZ"):
            p.inst(pyquil.gates.CZ(int(AC1[i]), int(AC2[i])))
    for i in range(0,nqubits):
        p.inst(pyquil.gates.MEASURE(i, ro[i]))
    p.wrap_in_numshots_loop(shots)
    return p

