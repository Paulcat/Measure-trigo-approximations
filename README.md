# Measure-trigo-approximations

This code illustrates the polynomial approximations of measures in the Wasserstein distance discussed in [[Catala, Hockmann, Kunis, Wageringel, 2022]](https://arxiv.org/abs/2203.10531v2). Given moments of a measure m up to some order n, we build polynomials of degree n such that the Wasserstein distance W(p,m) goes to zero as n increases. The distance is computed using a semidiscrete transport algorithm.

## Content
**examples:**
contains precomputed data for measures (either discrete, supported on a trigonometric curve or on a circle) and their moments.

**src:**
code for computing the polynomial approximations (proxys) and the semidiscrete transport (semiOT). The optimization is done using a bfgs implementation by M. Schmidt [[minFunc: unconstrained differentiable multivariate optimization in Matlab]](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)

## Usage
Simply run the example script ExampleProxy. Specified parameters can be tuned.
