# scourFoam

SedExnerFoam is an OpenFOAM solver in development aiming at the numerical study of scour around hydraulics structures like bridge piles. The objectiv is to simulate the hydrodynamics and both the bedload sediment transport and the suspension of sediment. Morphodynamics is handled with the Arbitrary Lagrangian Eulerian approach, which consists of a dynamic mesh getting distorted to track the morphology evolution.

- suspended load transport:

$$ \frac{\partial c_s}{\partial t} + \vec{\nabla}.[(\vec{u_s}+\vec{w_s})c_s] = \vec{\nabla}.(\epsilon_s\vec{\nabla} c_s) $$

- bedload transport and morphodynamics are described by the Exner equation:

$$ \frac{\partial z_b}{\partial t} + \vec{\nabla}.\vec{q_b}= D-E $$

The Exner equation is solved on a finite area mesh, which is a mesh made of finite polyhedral surfaces which can be curved in the 3D space.

# Repository

The repository is organized in different folders
- **solver**, all the files related to the solver
- **src**, all the librairies needed by the solver
- **tests**, cases to help identify bugs and missbehavior of model 
- **tutorials**, sedExnerFoam cases to test and validate the code

[OpenFoam](https://www.openfoam.com/) is the free, open source CFD software developed primarily by OpenCFD Ltd since 2004. The OpenFOAM trademark is owned by OpenCFD Ltd.
