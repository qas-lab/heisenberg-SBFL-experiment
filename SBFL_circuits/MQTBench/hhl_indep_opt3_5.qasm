// Benchmark created by MQT Bench on 2025-10-16
// For more info: https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: 2.0.1
// Qiskit version: 2.1.1
// Output format: qasm3

OPENQASM 3.0;
include "stdgates.inc";
gate qft_dg _gate_q_0, _gate_q_1, _gate_q_2 {
  swap _gate_q_0, _gate_q_2;
  h _gate_q_0;
  cp(-pi/2) _gate_q_1, _gate_q_0;
  h _gate_q_1;
  cp(-pi/4) _gate_q_2, _gate_q_0;
  cp(-pi/2) _gate_q_2, _gate_q_1;
  h _gate_q_2;
}
bit[1] c;
x $0;
h $1;
cp(11.366282220687872) $1, $0;
h $2;
cp(5.683141110343936) $2, $0;
h $3;
cp(2.841570555171968) $3, $0;
qft_dg $1, $2, $3;
cry(pi) $1, $4;
cry(pi) $2, $4;
cry(1.1714210053478946) $3, $4;
h $3;
cp(pi/2) $3, $2;
h $2;
cp(pi/4) $3, $1;
cp(pi/2) $2, $1;
h $1;
cp(-2.841570555171968) $1, $0;
h $1;
cp(-5.683141110343936) $2, $0;
h $2;
cp(-11.366282220687872) $3, $0;
h $3;
c[0] = measure $0;
