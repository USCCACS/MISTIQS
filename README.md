# MISTIQS: Multiplatform Software for Time-dependent Quantum Simulation
A full-stack, cross-platform software for generating, compiling, and executing quantum circuits for simulating quantum many-body dynamics of systems governed by Heisenberg Hamiltonians.

MISTIQS provides the following primary capabilities:

1) Generation of quantum circuits for performing quantum simulations of the dynamics of spins governed by input Heisenberg Hamiltonians

2) Translation of these quantum circuits into executable circuit objects for IBM, Rigetti, and Googele quantum devices (Qiskit circuit objects, Pyquil programs, and Cirq circuit objects respectively)

3) Compilation of circuits, either using the compilers native to the quantum computing platform of choice or using the domain-specific IBM and Rigetti compilers developed for TFIM simulations (more information here)

3) Execution of these circuits on IBM or Rigetti quantum processors

4) Post-processing and plotting of simulation results (limited to average magnetization data)

MISTIQS provides the user with extensive flexibility across its functionalities. Some examples include support for user-defined time dependence functions for external fields, 

To aid future researchers in benefiting from our domain-specific (DS) compilers, we have integrated them into the full stack quantum simulation software package provided in this repo.  This package provides the capability to i) generate circuits for the transverse-field Ising model (TFIM) simulations of interest in this article; ii) apply the DS compilers to these high-level circuits to generate hardware-executable circuits on either the Rigetti or IBM quantum processors; iii) execute the circuits on the Rigetti or IBM quantum processors; and iv) provide limited post-processing and plotting of results.  

If the user is solely interested in studying the performance of our DS compilers on their own TFIM circuits, they may also elect to bypass circuit generation and execution and have the software simply return the DS-compiled circuit.   The two provided tutorials provide demonstrations of how to carry out each use case within this software packge.
