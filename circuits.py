import cirq

q = cirq.LineQubit.range(8)

#circuit 1 - power law distribution
#uses independent Ry rotations encoding marginal probabilities (threshold=0.3)
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

#circuit 2 - Gaussian distribution
#uses uniformly controlled Ry gates w Mottonen encoding & threshold=0.07
circuit2 = cirq.Circuit([
    cirq.ry(1.576443177738)(q[0]),
    cirq.ry(1.574051774916)(q[1]),
    cirq.CX(q[0], q[1]),
    cirq.ry(0.445013839315)(q[1]),
    cirq.CX(q[0], q[1]),
    cirq.ry(1.572641510570)(q[2]),
    cirq.CX(q[0], q[2]),
    cirq.ry(0.240728633613)(q[2]),
    cirq.CX(q[0], q[2]),
    cirq.CX(q[1], q[2]),
    cirq.ry(0.118610301699)(q[2]),
    cirq.CX(q[1], q[2]),
    cirq.ry(1.571758138835)(q[3]),
    cirq.CX(q[0], q[3]),
    cirq.ry(0.123739135909)(q[3]),
    cirq.CX(q[0], q[3]),
    cirq.ry(1.571282726248)(q[4]),
    cirq.ry(1.571040234449)(q[5]),
    cirq.ry(1.570918369802)(q[6]),
    cirq.ry(1.570857359468)(q[7]),
    cirq.measure(*q, key='m')
])

#circuit 3 - linearly decreasing distribution
# uses uniformly controlled Ry gates w Mottonen encoding & threshold=0.1
circuit3 = cirq.Circuit([
    cirq.ry(1.147579938763)(q[0]),
    cirq.ry(1.319792002067)(q[1]),
    cirq.CX(q[0], q[1]),
    cirq.ry(0.104920473487)(q[1]),
    cirq.CX(q[0], q[1]),
    cirq.ry(1.433512853812)(q[2]),
    cirq.ry(1.499501033388)(q[3]),
    cirq.ry(1.534690139035)(q[4]),
    cirq.ry(1.552677501555)(q[5]),
    cirq.ry(1.561728320221)(q[6]),
    cirq.ry(1.566261235792)(q[7]),
    cirq.measure(*q, key='m')
])