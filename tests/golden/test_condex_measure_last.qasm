version 1.0
# this file has been automatically generated by the OpenQL compiler please do not modify it manually.
qubits 7

.kernel_measure
    measure q[0], b[0]
    wait 14
    { cond(b[0]) measure q[1], b[1] | cond(!b[0]) measure q[2], b[2] }
    wait 14
    { cond(b[1]) x q[3] | cond(b[2]) y q[5] }
