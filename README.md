# Controlling Epidemics via Mean Field Games

## About

This github repository contains the C++ and python codes for the epidemics mean-field games algorithm from ["Controlling Propagation of epidemics via mean-field games" ](https://arxiv.org/abs/2006.01249). 

## Prerequisite

1. C++ compiler (the codes use C++14)
2. FFTW library. ([link](http://www.fftw.org/fftw2_doc/fftw_6.html) to the instruction for the installation)
3. python3 with numpy and matplotlib

## How to run

1. Run the following command in the terminal to compile the codes:

```make```

2. After compilation is done, run the following command to run the algorithm:

```./SIR-PDHG 32 32 12 0.1 0.1 5e-7 1000 50 0.6 0.2```

3. Run the following python codes to plot the graph and create the video of the solutions from the algorithm.

```python plot-contour.py```
