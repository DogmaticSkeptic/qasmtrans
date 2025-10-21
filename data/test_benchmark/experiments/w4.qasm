OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];

// Create |W4> = (|1000>+|0100>+|0010>+|0001>)/2

ry(1.0471975512) q[0];
cx q[0], q[1];
ry(1.2309594173) q[1];
cx q[1], q[2];
ry(1.5707963268) q[2];
cx q[2], q[3];
