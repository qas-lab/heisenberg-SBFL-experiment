// Benchmark created by MQT Bench on 2025-10-16
// For more info: https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 2.0.1
// Qiskit version: 2.1.1
// Output format: qasm3

OPENQASM 3.0;
include "stdgates.inc";
bit[5] meas;
qubit[5] q;
cx q[1], q[3];
cx q[1], q[4];
ccx q[0], q[2], q[1];
ccx q[1], q[3], q[4];
cx q[1], q[3];
ccx q[0], q[2], q[1];
cx q[0], q[2];
cx q[1], q[3];
barrier q[0], q[1], q[2], q[3], q[4];
meas[0] = measure q[0];
meas[1] = measure q[1];
meas[2] = measure q[2];
meas[3] = measure q[3];
meas[4] = measure q[4];
