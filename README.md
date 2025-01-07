# SedExnerFoam

SedExnerFoam is an OpenFOAM solver in development aiming at the numerical study of various sediment transport and morphodynamics problems. It is based on the solver pimpleFoam for te hydrodynamics. The sediments can be either transported in suspension in the fluid or move but stay in contact with the bed. Morphodynamics is handled with the Arbitrary Lagrangian Eulerian approach. With this approach, computational domain stops at the sediment bed and thus does not include the granular material below. The sediments movements along this boundary (bedload) or through it (erosion/deposition) drive its morphology evolution. The mesh is dynamic and get distorted in time so that the bed boundary match the geometry changes.

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
- **tutorials**, examples of sedExnerFoam cases and applications with validation

[OpenFoam](https://www.openfoam.com/) is the free, open source CFD software developed primarily by OpenCFD Ltd since 2004. The OpenFOAM trademark is owned by OpenCFD Ltd.

