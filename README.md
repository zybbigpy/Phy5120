# Phy5120

## How to build my code

First of all, clone my code from Github.

```bash
git clone git@github.com:zybbigpy/Phy5120.git
```

Then, use `git pull` to get the latest code and report.

```bash
git pull origin master
```

If your system support `Make`, just use the command below. if you have `pandoc` and LaTeX in your system, you can build the pdf file. Otherwise you just enter the report directory to get the pdf file. If you want to change the parameter, read the `Makefile` and you will find it is so easy.

```bash
# go to hw1
cd hw1
# run code
make run
# make pdf file
make doc
# clean data
make clean
# format code
make format
```
