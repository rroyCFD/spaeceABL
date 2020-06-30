{
    // kappat = turbulence->nut()/Prt;
    // kappat.correctBoundaryConditions();

    // volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + kappat);
    volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + turbulence->nut()/Prt);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
      - fvc::div(qwall)
      - temperatureSourceTerm.force()
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    Boussinesq.updateDensityField();
}
