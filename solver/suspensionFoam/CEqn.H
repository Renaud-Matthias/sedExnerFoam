{
    volScalarField epsEff("epsEff", turbulence->nu());
    
    fvScalarMatrix CEqn
    (
          fvm::ddt(C)
        + fvm::div(phi, C)
        - fvm::laplacian(epsEff, C)
        ==
          fvOptions(C)
    );
    
    CEqn.relax();
    
    fvOptions.constrain(CEqn);
    
    CEqn.solve();
    
    fvOptions.correct(C);
}
