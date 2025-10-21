OPENQASM 2.0;
include "qelib1.inc";

qreg q[4];
creg c[4];

// Layer 1
ry(0.442699) q[0];
rz(0.563599) q[0];
ry(0.462699) q[1];
rz(0.578599) q[1];
ry(0.482699) q[2];
rz(0.593599) q[2];
ry(0.502699) q[3];
rz(0.608599) q[3];
// Entangling layer 1
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 2
ry(0.492699) q[0];
rz(0.603599) q[0];
ry(0.512699) q[1];
rz(0.618599) q[1];
ry(0.532699) q[2];
rz(0.633599) q[2];
ry(0.552699) q[3];
rz(0.648599) q[3];
// Entangling layer 2
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 3
ry(0.542699) q[0];
rz(0.643599) q[0];
ry(0.562699) q[1];
rz(0.658599) q[1];
ry(0.582699) q[2];
rz(0.673599) q[2];
ry(0.602699) q[3];
rz(0.688599) q[3];
// Entangling layer 3
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 4
ry(0.592699) q[0];
rz(0.683599) q[0];
ry(0.612699) q[1];
rz(0.698599) q[1];
ry(0.632699) q[2];
rz(0.713599) q[2];
ry(0.652699) q[3];
rz(0.728599) q[3];
// Entangling layer 4
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 5
ry(0.642699) q[0];
rz(0.723599) q[0];
ry(0.662699) q[1];
rz(0.738599) q[1];
ry(0.682699) q[2];
rz(0.753599) q[2];
ry(0.702699) q[3];
rz(0.768599) q[3];
// Entangling layer 5
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 6
ry(0.692699) q[0];
rz(0.763599) q[0];
ry(0.712699) q[1];
rz(0.778599) q[1];
ry(0.732699) q[2];
rz(0.793599) q[2];
ry(0.752699) q[3];
rz(0.808599) q[3];
// Entangling layer 6
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 7
ry(0.742699) q[0];
rz(0.803599) q[0];
ry(0.762699) q[1];
rz(0.818599) q[1];
ry(0.782699) q[2];
rz(0.833599) q[2];
ry(0.802699) q[3];
rz(0.848599) q[3];
// Entangling layer 7
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 8
ry(0.792699) q[0];
rz(0.843599) q[0];
ry(0.812699) q[1];
rz(0.858599) q[1];
ry(0.832699) q[2];
rz(0.873599) q[2];
ry(0.852699) q[3];
rz(0.888599) q[3];
// Entangling layer 8
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 9
ry(0.842699) q[0];
rz(0.883599) q[0];
ry(0.862699) q[1];
rz(0.898599) q[1];
ry(0.882699) q[2];
rz(0.913599) q[2];
ry(0.902699) q[3];
rz(0.928599) q[3];
// Entangling layer 9
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Layer 10
ry(0.892699) q[0];
rz(0.923599) q[0];
ry(0.912699) q[1];
rz(0.938599) q[1];
ry(0.932699) q[2];
rz(0.953599) q[2];
ry(0.952699) q[3];
rz(0.968599) q[3];
// Entangling layer 10
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];
cx q[3],q[0];

// Measurement
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
