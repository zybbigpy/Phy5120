import numpy as np
import math
from matplotlib import pyplot as plt


length = 8
g2m_bins = 33
m2k_bins = 33
k2g_bins = 33

def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n


line = ''
k1 = np.linspace(0, 0.5, g2m_bins)
# print(k)
print("set wf_dyn JD # JD may be better for unoccupied states\nset nempty 8 # add 8 empty bands")

line+="set wf_dyn JD # JD may be better for unoccupied states\nset nempty 8 # add 8 empty bands\n"
# line+="kpoint delete 0 0 0\n"

print("#","-"*30,"Gamma point to M point","-"*30)
line+="# ------------------------------ Gamma point to M point ------------------------------\n"
for i in k1:
    print("kpoint add%11.8f%11.8f%11.8f%8.3f" % (truncate(i,length), 0, 0, 0))
    line+="kpoint add%11.8f%11.8f%11.8f%8.3f\n" % (truncate(i,length), 0, 0, 0)

print("#","-"*30,"M point to K point","-"*30)
line+="# ------------------------------ M point to K point ------------------------------\n"
k2 = np.linspace(0.5, 1.0/3, m2k_bins)
k3 = np.linspace(0, 1.0/3, m2k_bins)

for i in range(len(k2)):
    print("kpoint add%11.8f%11.8f%11.8f%8.3f" % (truncate(k2[i],length), truncate(k3[i],length), 0, 0))
    line+="kpoint add%11.8f%11.8f%11.8f%8.3f\n" % (truncate(k2[i],length), truncate(k3[i],length), 0, 0)


print("#","-"*30,"K point to Gamma point","-"*30)
line+="# ------------------------------ K point to Gamma point ------------------------------\n"
k4 = np.linspace(1.0/3, 0, k2g_bins)
k5 = np.linspace(1.0/3, 0, k2g_bins)
for i in range(len(k4)):
    print("kpoint add%11.8f%11.8f%11.8f%8.3f" % (truncate(k4[i],length), truncate(k5[i],length), 0, 0))
    line+="kpoint add%11.8f%11.8f%11.8f%8.3f\n" % (truncate(k4[i],length), truncate(k5[i],length), 0, 0)


print("#","-"*50)
line+="# -------------------------------------------------------------------------------------\n"
print("run 0 1 200 # 200 non-SCF steps, which may not be enough!\nrun 0 1 # double check the convergence of eigenvalues")

line+="run 0 1 300 # 200 non-SCF steps, which may not be enough!\nrun 0 1\nsave prob2.xml\n# double check the convergence of eigenvalues\n#\n"

with open('band.i','w') as f:
    f.write(line)