OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];

// superposition on first two qubits
h q[0];
h q[1];

// placeholder for controlled modular multiplication (using simple CX chain)
cx q[0], q[2];
cx q[0], q[3];
cx q[1], q[2];
cx q[1], q[3];

// inverse QFT on control register (q0,q1)
cu1(-pi/2) q[0], q[1];
h q[1];
h q[0];
