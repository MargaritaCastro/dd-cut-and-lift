# About This Work

This repository includes the code and data instances used in the paper: 

*Castro, M.P., Cire, A.A. & Beck, J.C. A combinatorial cut-and-lift procedure with an application to 0–1 second-order conic programming. Math. Program. (2021).*

A preliminary version of this paper can be found in [arXiv](https://arxiv.org/abs/2003.06363v2). 

A final version of the paper is availabale [here](https://doi.org/10.1007/s10107-021-01699-y).
 

# The Code

The main portion of the code (i.e., BDD cuts and lifting) is written in C++ (std 14). 
We also provide the data generator, which is written in  Python 3


## Files includes

* makefile
* data instances used in the paper: SOC knapsack and SOC General
* data instances generator
* C++ source code 


## Code requirements

* CPLEX 12.9 -- the code probably works for newer versions of CPLEX but we haven't test it
* Boost library fdor C++ [(here)](https://www.boost.org/)
* C++ with std=14
* Python 3 (only for instance generation)


## Running instructions

Please refer to the main file (i.e., src/main.cpp) to see the code usage and meaning of flags. 

## Support
Feel free to email me at `margarita.castro[at]ing.puc.cl` if you have any questions or find any issues. 

