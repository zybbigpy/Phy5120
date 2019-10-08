# Homework 2, Miao Wangqian

## Molecular dynamics simulations of Lennard-Jones Argon



```bash
# use the following command in bash
# you will get all figures and .gro file
# in the directory, the value of D shown
# on your screen
python md.py
```

## 1

The temperature in the reduced units is $\frac{203K}{119.8K}=1.70$

The number density in reduced units is $\frac{4}{1(\sigma)} =4$

## 2

The total energy is conserved, but kinetic energy and potential energy will fluctuate.

![energy](../figs/energy_time.png)

## 3

![md](../Argon_MD.gif)

If the time step is very small, the change of position value will be small. so we do not need to write down every time step.

## 4

The radial distribution function is shown:

![radial](../figs/raditial_func.png)

## 5

It should be the same when we use the two different methods. But my numerical results show that it has some difference. When computing MSDs, we should use unfolded coordinates because we want to get the true displacement of the atoms.

There may be something wrong in my code ?

```python
# in every production step,
# I use the following code to record msd
msd = np.sum((pos - pos_0)**2) / N
msd_list.append(msd)
# calculate the gradient of msd_list
D1 = np.average(np.gradient(np.array(msd_list))/ (0.004*6)
```

