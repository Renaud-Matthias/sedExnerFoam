// Volume-to surface mapping object
const volSurfaceMapping vsm(aMesh);

volScalarField Zb
(
    IOobject
    (
        "Zb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, Zero)
);
