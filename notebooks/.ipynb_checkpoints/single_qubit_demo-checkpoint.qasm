OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];

// Single-qubit-only program for 8 qubits (no two-qubit gates)
// Sequence: X then RZ(pi/2) on each qubit, then measure.
x q[0];
rz(1.57079632679) q[0];

x q[1];
rz(1.57079632679) q[1];

x q[2];
rz(1.57079632679) q[2];

x q[3];
rz(1.57079632679) q[3];

x q[4];
rz(1.57079632679) q[4];

x q[5];
rz(1.57079632679) q[5];

x q[6];
rz(1.57079632679) q[6];

x q[7];
rz(1.57079632679) q[7];

measure q -> c;
