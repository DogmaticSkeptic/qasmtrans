OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];
creg c[4];

// Hardware-efficient ansatz layer
ry(0.4) q[0];
ry(0.6) q[1];
ry(0.5) q[2];
ry(0.7) q[3];

rz(0.1) q[0];
rz(0.2) q[1];
rz(0.3) q[2];
rz(0.4) q[3];

cx q[0], q[1];
cx q[1], q[2];
cx q[2], q[3];

ry(0.2) q[0];
ry(0.3) q[1];
ry(0.4) q[2];
ry(0.5) q[3];

cx q[2], q[1];
cx q[3], q[2];
cx q[1], q[0];

measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
