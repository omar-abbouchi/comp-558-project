from circuits import circuit1, circuit2, circuit3
from score_circuits import analyze_circuits

score = analyze_circuits(circuit1, circuit2, circuit3)
print('Score:', score)