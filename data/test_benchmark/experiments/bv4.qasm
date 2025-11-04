OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];
creg c[4];

// initialise ancilla in |->
x q[3];
h q[3];

// input qubits in superposition
h q[0];
h q[1];
h q[2];

// oracle for string 1011 (bits on q0,q1,q2)
cx q[0], q[3];
cx q[1], q[3];
cx q[2], q[3];

// uncompute ancilla
h q[3];

// Hadamard to read the string
h q[0];
h q[1];
h q[2];

measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
