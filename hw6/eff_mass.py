from lxml import etree
import re
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


plt.rc('text', usetex=True)
# Computer Modern font
plt.rc('font', family='serif')

file = open("out.o", 'rb')
content = file.read()
file.close()
xml = etree.XML(content)
k_eigen_values = xml.xpath("//eigenvalues")
k_eigen_values = [[float(j) for j in i.text.split()] for i in k_eigen_values][0:96]

# eng[13] is the bottom of the conduction band
y = [eng[13] for eng in k_eigen_values]
# from M to K path eV to J
y = np.array(y[33:66])*1.6e-19


k1 = np.linspace(0.5, 1.0/3, 33)
k2 = np.linspace(0, 1.0/3, 33)

path = [np.array([k1[i], k2[i], 0]) - np.array([0.5, 0, 0]) for i in range(33)]

val = [k[0]*np.array([1.041349, 0.601226, -0.000000]) + k[1]*
np.array([0.000000, 1.202431, 0.000000])  for k in path]

val_abs = [np.sqrt(i[0]**2+i[1]**2)/0.529e-10 for i in path]

print(len(y))
print(len(val_abs))

def func(x, a, b, c):

    return a*x*x+b*x+c

h_bar = 6.62607015/np.pi/2*1e-34
popt, pcov = curve_fit(func, val_abs, y)
a, b, c = popt[0], popt[1], popt[2]
mass = h_bar**2/(2*a)
print("Effective mass is: ", mass)
yvals = func(np.array(val_abs), *popt)
y1 = np.array(y)/1.6e-19
yvals = np.array(yvals)/1.6e-19
fig, ax = plt.subplots()
ax.plot(val_abs, y1, '.', label='original values')
ax.plot(val_abs, yvals, 'r', label='fitting values')
ax.set_xticks([0,val_abs[-1]])
ax.set_xlim(0, val_abs[-1])
ax.set_xticklabels(["$M$", "K"])
ax.axvline(x=0, color="black")
ax.axvline(x=val_abs[-1], color="black")
ax.set_ylabel("Energy(eV)")
# plt.xticks_labels("$\Gamma$","M")
ax.legend()
plt.savefig("./effmass.png", dpi =500)