import numpy as np 
import matplotlib.pyplot as plt


x = np.array(np.arange(20, 140, 20))
y = np.array([-16.08873745, -16.92295762, -17.11900488, -17.15675309, -17.16309513, -17.16403754])
print(x, y)
plt.xlabel("ecut:(Ry)")
plt.ylabel("total energy (Hartree)")
plt.plot(x, y)
plt.savefig("./1-b.png")