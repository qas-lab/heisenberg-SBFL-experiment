// Benchmark created by MQT Bench on 2025-10-16
// For more info: https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 2.0.1
// Qiskit version: 2.1.1
// Output format: qasm3

OPENQASM 3.0;
include "stdgates.inc";
bit[5] meas;
h $4;
cp(pi/2) $4, $3;
h $3;
cp(pi/4) $4, $2;
cp(pi/2) $3, $2;
h $2;
cp(pi/8) $4, $1;
cp(pi/4) $3, $1;
cp(pi/2) $2, $1;
h $1;
cp(pi/16) $4, $0;
cp(pi/8) $3, $0;
cp(pi/4) $2, $0;
cp(pi/2) $1, $0;
h $0;
barrier $4, $3, $2, $1, $0;
meas[0] = measure $4;
meas[1] = measure $3;
meas[2] = measure $2;
meas[3] = measure $1;
meas[4] = measure $0;
