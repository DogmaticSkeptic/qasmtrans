OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];

// Treat q0,q1 as input number, q2 as carry, q3 as target
// Initialize to |q0 q1 q2 q3> = |00 01> (for demonstration)
x q[0];

// ripple-carry increment by 1
cx q[0], q[1];
ccx q[0], q[1], q[2];
cx q[0], q[1];
cx q[2], q[3];
