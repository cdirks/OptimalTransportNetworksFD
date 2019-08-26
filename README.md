# Optimal transportation networks via Mumford-Shah-type image inpainting 

## General information 

The code implements the method proposed in [1](), Chapter 4.3, for solving the branched transport [2](http://www.uvm.edu/pdodds/research/papers/others/2003/xia2003a.pdf),[3](https://pdfs.semanticscholar.org/d766/7ac83e8dd7c8ce452fe63775a3ddd705efd9.pdf) and urban planning [4](http://www.numdam.org/article/COCV_2005__11_1_88_0.pdf) problem via a convex reformulation as a Mumford–Shah-type image inpainting problem introduced by [5](https://arxiv.org/abs/1601.07402).

The implementation is based on [Matlab](https://www.mathworks.com/products/matlab.html). 


## Run example 

To run the code, open Matlab and execute Example.m. This will compute a simple example of branched transport from one point source to two point sinks with equal mass. For an other example or a different setting, the following parameters can be changed within the example file:

	type: Example type (BT or UP)
	n: Image size (n by n grid) 
	savename: Savename 
	maxIter: Maximal number of primal-dual iterations
	tol: Primal-dual iteration error tolerance 
	maxIterProj: Maximal number of Dykstra projection iterations 
	tolProj: Error tolerance for Dykstra iteration 
	alpha: Branched transport parameter 
	a: Urban planning parameter
	epsilon: Urban planning parameter 
	example: Example (see getExample.m)


[1] Carolin Dirks. Numerical methods for transportation networks. PhD thesis, 2019.

[2] Qinglan Xia. Optimal paths related to transport problems. *Commun. Contemp. Math.*, 5(2):251--279, 2003. 

[3] Francesco Maddalena, Sergio Solimini, and Jean-Michel Morel. A variational model of irrigation patterns. *Interfaces Free Bound.*, 5(4):391--415, 2003.

[4] Alessio Brancolini and Giuseppe Buttazzo. Optimal networks for mass transportation problems. *ESAIM Control Optim. Calc. Var.*, 11(1):88--101 (electronic), 2005.

[5] Alessio Brancolini, Carolin Rossmanith, Benedikt Wirth. Optimal micropatterns in 2D transport networks and their relation to image inpainting. *Archive for Rational Mechanics and Analysis*, 228(1):279--308, 2018.

