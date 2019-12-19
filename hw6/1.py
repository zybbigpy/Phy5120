import numpy as np
from matplotlib import pyplot as plt
k = np.array([1,2,3,4,5])

e = np.array([    -88.94624379,
    -89.24607856,
    -89.21231193,
    -89.06546033,
    -89.23481561])

plt.title("Total Energy with different kpoint mesh(n x n x 1)")
plt.ylabel("Energy (Hartree)")
plt.xlabel("n(grid number)")
plt.plot(k,e)
plt.savefig("./1.png",dpi=400)