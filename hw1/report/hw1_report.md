# PHYS ​5120: Homework 1

Author: MIAO Wangqian, Student ID: 20617902

## 1. How to build my code

First of all, clone my code from Github.

```bash
git clone git@github.com:zybbigpy/Phy5120.git
```

Then, use `git pull` to get the latest code and report.

```bash
git pull origin master
```

If your system support `Make`, just use the command below. if you have `pandoc` and LaTeX in your system, you can build the pdf file. Otherwise you just enter the report directory to get the pdf file.

```bash
# go to hw1
cd hw1
# run code
make run
# make pdf file
make doc
# clean data
make clean
```

## 2. The Linear and nonlinear pendulums

### 2.1 Solution

The equation of motion is in the format of:

$$ \frac{\mathrm{d^2}\theta}{\mathrm{d} t^2} + \frac{g}{\ell}\theta = 0$$

And the solution of the differential equation is:

$$\theta =  A \cos (\sqrt\frac{g}{\ell} t + \delta)$$

There are two parameters $A, \delta$ in the solution because we do not know the initial condition $\theta(t=0), \dot{\theta}(t=0)$.

The swing period is:

$$T = 2\pi \sqrt\frac{\ell}{g}$$

### 2.2 Solution




