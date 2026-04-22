import cirq

q = cirq.LineQubit.range(8)

#circuit 1 power law distribution
#uses independent Ry rotations encoding marginal probabilities & threshold=0.3
circuit1 = cirq.Circuit([
    cirq.ry(1.699170175247)(q[0]),
    cirq.ry(1.651263732269)(q[1]),
    cirq.ry(1.619029558502)(q[2]),
    cirq.ry(1.598695347295)(q[3]),
    cirq.ry(1.586447693226)(q[4]),
    cirq.ry(1.579330085044)(q[5]),
    cirq.ry(1.575323294140)(q[6]),
    cirq.ry(1.573139094627)(q[7]),
    cirq.measure(*q, key='m')
])

#circuit 2 gaussian distribution optimized using Nelder-Mead
# hardware-efficient variational ansatz: Ry layer -> CNOTs on MSB -> Ry layer

circuit2 = cirq.Circuit([
    cirq.ry(1.572976433791605)(q[0]),
    cirq.ry(3.1405035674334583)(q[1]),
    cirq.ry(1.7758615236362747)(q[2]),
    cirq.ry(1.4745930518771186)(q[3]),
    cirq.ry(7.2951958147507945)(q[4]),
    cirq.ry(0.2712157238090181)(q[5]),
    cirq.ry(2.9307086288795716)(q[6]),
    cirq.ry(0.21004458915148114)(q[7]),
    cirq.CX(q[0], q[1]),
    cirq.CX(q[0], q[2]),
    cirq.CX(q[0], q[3]),
    cirq.ry(3.698500532273342)(q[0]),
    cirq.ry(1.4728888016414534)(q[1]),
    cirq.ry(3.1465217118268507)(q[2]),
    cirq.ry(6.289320266437411)(q[3]),
    cirq.ry(0.556806901392646)(q[4]),
    cirq.ry(1.2963811791653195)(q[5]),
    cirq.ry(1.7813127306250236)(q[6]),
    cirq.ry(4.5022871177937525)(q[7]),
    cirq.measure(*q, key='m')
])

#circuit 3 linearly decreasing distribution
# a hardware-efficient variational ansatz: Ry layer -> CX(q0,q1) -> Ry layer
#also optimized using Nelder-Mead on exact statevector TVD
circuit3 = cirq.Circuit([
    cirq.ry(0.6521204370732243)(q[0]),
    cirq.ry(4.872169979438421)(q[1]),
    cirq.ry(2.48889833908835)(q[2]),
    cirq.ry(0.06919625454162216)(q[3]),
    cirq.ry(4.416622127788479)(q[4]),
    cirq.ry(3.819420563447241)(q[5]),
    cirq.ry(0.3527756995560756)(q[6]),
    cirq.ry(2.5683015058585426)(q[7]),
    cirq.CX(q[0], q[1]),
    cirq.ry(5.780429780674856)(q[0]),
    cirq.ry(2.8079941569341855)(q[1]),
    cirq.ry(5.281631723744741)(q[2]),
    cirq.ry(4.683526054646446)(q[3]),
    cirq.ry(0.3148016107345143)(q[4]),
    cirq.ry(4.024946183050365)(q[5]),
    cirq.ry(1.2131661352707823)(q[6]),
    cirq.ry(2.1465240588181382)(q[7]),
    cirq.measure(*q, key='m')
])