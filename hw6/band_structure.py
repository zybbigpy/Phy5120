from lxml import etree
import re
import numpy as np
from matplotlib import pyplot as plt
g2m_bins = 33
m2k_bins = 33
k2g_bins = 33

plt.rc('text', usetex=True)
# Computer Modern font
plt.rc('font', family='serif')
file = open("out.o", 'rb')
content = file.read()
file.close()
xml = etree.XML(content)
k_eigen_values = xml.xpath("//eigenvalues")
k_eigen_values = [[float(j) for j in i.text.split()] for i in k_eigen_values][0:96]
# for i in k_eigen_values:
#     print(i)
x = np.arange(len(k_eigen_values))
print("Band gap at Gamma point is: %.5feV" % (k_eigen_values[64][13]-k_eigen_values[64][12]))
fig, ax = plt.subplots()

for i in range(12, 14):
    y = []
    for j in k_eigen_values:
        y.append(j[i])
    plt.plot(x, y)
ax.set_xticks([0,g2m_bins-1,g2m_bins+m2k_bins-2,g2m_bins+m2k_bins+k2g_bins-4])
ax.set_xticklabels(["$\Gamma$","M","K","$\Gamma$"])
ax.set_xlim(0, max(x))
ax.set_ylabel("Engergy (eV)")
ax.set_title("Band Structure of MoS$_2$")
ax.axvline(x=0, color="black")
ax.axvline(x=g2m_bins-1, color="black")
ax.axvline(x=g2m_bins+m2k_bins-2, color="black")
ax.axvline(x=g2m_bins+m2k_bins+k2g_bins-4, color="black")
plt.show()
#plt.savefig("./100.png", dpi =500)

# print(content)
