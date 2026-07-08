// Benchmark created by MQT Bench on 2025-10-16
// For more info: https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 2.0.1
// Qiskit version: 2.1.1
// Output format: qasm3

OPENQASM 3.0;
include "stdgates.inc";
bit[4] c;
h $0;
h $1;
h $2;
h $3;
x $4;
cp(pi/16) $4, $0;
cp(pi/8) $4, $1;
cp(pi/4) $4, $2;
cp(pi/2) $4, $3;
h $3;
cp(-pi/2) $2, $3;
cp(-pi/4) $1, $3;
cp(-pi/8) $0, $3;
h $2;
cp(-pi/2) $1, $2;
cp(-pi/4) $0, $2;
h $1;
cp(-pi/2) $0, $1;
h $0;
barrier $3, $2, $1, $0, $4;
c[0] = measure $3;
c[1] = measure $2;
c[2] = measure $1;
c[3] = measure $0;
