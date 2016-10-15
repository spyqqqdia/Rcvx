## Rcvx ##

The goal of this project is to implement some interior point solvers 
for dense convex optimization.
The algorithms follow the description in Boyd and Vandenberghe, 
[Convex Optimzation](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf).

Here we are prototyping what is intended to become a Scala implementation later.
For savvy R programmers there is little of interest here.
For beginners the following might be of some interest:

+ Use of S4 classes, derivation in R.
+ All manner of linear algebra manipulations, see [nullspace](docs/nullspace.pdf)
+ Example of how a multifile project might be organized in R.

With current documentation the algorithms for convex optimization will be
incomprehensible. It is intended to rectify this later.

Why do we need an implementation of these algorithms in Scala when for example  
the [real experts](https://github.com/cvxgrp/) have treated the problem in   
[greater generality](https://github.com/cvxgrp/scs) and provide a scala interface?  

One answer is that in order to use these tools you need to know what the _partial order_  
induced by a convex _cone_ is and how, using a product of such cones, you can formulate   
all manner of interesting convex optimization problems. 

Here we are dealing only with optimization problems of the form

    minimize f(x) subject to $g_i(x)\leq 0$ and Ax=b

with twice continuously differentiable, _convex_ functions $f,g_i:R^n\to R$ and the inequality
$g_i(x)\leq 0$ uses only the order on the real numbers.

It should be easy to set up such a problem, but in practice it is not. The algorithms all
work on arbitrary convex functions subject to the above restrictions but existing R-packages
or other software often limit the objective function f(x) and constraint functions g_i(x)
to be of simple form, e.g. linear or quadratic or maybe defined in some modelling language.

Even though this latter approach can be fairly general it will not encompass functions 
the value of which is the result of running a program for example. It is intended to make this 
possible and reasonably easy by means of virtual functions and class derivation.

No attempt will be made to deal with sparse matrices in a special way.
This means that the problem size must be small (say up to 2000 decision variables
and similar numbers of constraints). There are enough interesting problem which are
in the scope of such an aproach.


 

