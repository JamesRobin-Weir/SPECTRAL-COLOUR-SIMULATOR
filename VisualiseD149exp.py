import matplotlib.pyplot as plt
import numpy as np

file_path="D149_exp"

Lambda=[]
RSP=[]

with open(file_path) as file:
    for _ in range(1):
        next(file)
    for line in file:
        a, b = line.rstrip().split(" ", 1)
        Lambda.append(float(a))
        RSP.append(float(b))


plt.plot(Lambda, RSP)
plt.title("Spectre de transmission du D149 experimental")
plt.xlabel("Longeur d'onde (nm)")
plt.ylabel('Transmittance du matériaux')
plt.show()
