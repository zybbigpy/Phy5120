import numpy as np
import math

# three primitive vetors in the space (unit A  1 A = 1.8897259886 Bohr)
convert = 1.8897259886
a = []

# read in the position in file
pos_file = open("./MoS2pos.txt")
atomname = []
xyz = []
newxyz = []

for line in pos_file:
    pos = line.split()
    index = 0
    if(len(pos) == 4):
        xyz.append(np.array( [float(i) for i in pos[0:3]]))
        atomname.append(pos[3])
    elif(len(pos) == 3):
        a.append(np.array([float(i) for i in pos])*convert)
        

for vec in a:
    #print(vec)
    print("%10.7f %10.7f %10.7f"%(vec[0], vec[1],vec[2]))

print("+"*30)

for atoms in xyz:
    #print(atoms)
    #print("-"*30)
    #print("%10.7f %10.7f %10.7f"%(atoms[0], atoms[1],atoms[2]))
    newpos = np.array(atoms[0]*a[0]+atoms[1]*a[1]+atoms[2]*a[2])
    print("%10.7f %10.7f %10.7f"%(newpos[0], newpos[1],newpos[2]))
    #print("-"*30)
#print(xyz)
print(atomname)
#print(newxyz)
pos_file.close()