import matplotlib.pyplot as plt
import numpy as np

file_path="/home/arthur/Documents/FL 03 – Simulation numérique de la couleur d’un colorant-20250908/oeil"

Lambda=[]
x_bar=[]
y_bar=[]
z_bar=[]

with open(file_path) as file:
    for _ in range(1):
        next(file)
    for line in file:
        a, b, c, d = line.rstrip().split("	")
        Lambda.append(float(a))
        x_bar.append(float(b))
        y_bar.append(float(c))
        z_bar.append(float(d))


plt.plot(Lambda, x_bar, label="x_bar")
plt.plot(Lambda, y_bar, label="y_bar")
plt.plot(Lambda, z_bar,label="z_bar")

plt.title("Spectre de reponse de l'oeil")
plt.xlabel("Longeur d'onde (nm)")
plt.ylabel("Reponse de l'oeuil")
plt.legend(loc="upper right")
plt.show()
