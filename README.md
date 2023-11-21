# scourFoam

ScourFoam is an OpenFOAM solver in development aiming at the numerical study of scour around hydraulics structures like bridge piles. The objectiv is to simulate the hydrodynamics and both the bedload sediment transport and the suspension of sediment. Morphodynamics is handled with the Arbitrary Lagrangian Eulerian approach, which consists of a dynamic mesh getting distorted to follow the morphology evolution.
- suspended load transport:
$$
\frac{\partial c}{\partial t} + \vec{\nabla}.[(\vec{u}+\vec{w_s})c] = \vec{\nabla}.(\epsilon_s\vec{\nabla} c)
$$
- bedload transport and morphodynamics are described by Exner equation:
$$
\frac{\partial z_b}{\partial t} + \vec{\nabla}.\vec{q_b}= D-E
$$
The Exner equation is solved on a finite area mesh, which is a mesh made of finite polyhedral surfaces which can be curved in the 3D space.

# Repository

The repository is organized in different folders
- tutorials - various OpenFoam cases to test and validate the code
- solver - all the files related to the solver
- src - all the librairies needed by the solver

[OpenFoam](https://www.openfoam.com/) is the free, open source CFD software developed primarily by OpenCFD Ltd since 2004. The OpenFOAM trademark is owned by OpenCFD Ltd.
