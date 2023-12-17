# C++ library of the paper: “Early Exercise Decision in American Options with Dividends, Stochastic Volatility and Jumps ”.

Code comes with no warranty. If you use parts of it, please cite the original source.

The library relies on other available open source libraries:
You will need
libnewmat : http://www.robertnz.net/nm_intro.htm
it’s a simple matrix algebra matrix
libfftw3.a: for fast fourier transform.  http://www.fftw.org
librecipes.a : not sure it’s strictly necessary. It’s a library of c++ codes for finance. You should only need the header that we included in the 'include' folder. Here is the source http://finance.bi.no/~bernt/gcc_prog/recipes/

As an example, we provide the main_pricing.cpp. The code calls the libFRP library and computes the American option prices under the Heston, Merton and Black Scholes models. The user should change the parameters in the box  /*INPUT PARAMETERS*/ at the beginning of the main.

The code has been tested on both Windows and Mac. If you compile the code from the command line in a Linux-like environment (tested on mac) you should do the following steps

_Compile the library libFRP.so_

make lib

_Compile the main_pricing.cpp example, and link it to the library_

make exec

_Run the executable main_pricing_

make run
