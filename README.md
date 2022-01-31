# DDA.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://fekad.github.io/DDA.jl/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://fekad.github.io/DDA.jl/dev)
[![Run tests](https://github.com/fekad/DDA.jl/actions/workflows/test.yml/badge.svg)](https://github.com/fekad/DDA.jl/actions/workflows/test.yml)


Discrete Dipole Approximation


## Contributions

- Adam Fekete (adam.fekete@unamur.be)

This work was carried out at UNamur in Luc Henrard's group . The authors acknowledge financial support from the Communauté Française de Belgique through the Action de Recherche Concertée (ARC - No. 19/24-102)
SURFASCOPE - Surface Enhanced Spectroscopy by Second-Principles Calculations.



solve(prob::ODEProblem,alg;kwargs)

Solves the ODE defined by prob using the algorithm alg. If no algorithm is given, a default algorithm will be chosen.

Solution Handling