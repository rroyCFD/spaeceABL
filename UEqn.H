    // Solve the momentum equation
    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(U) 
      + fvc::div(phi, U)
      + MRF.DDt(U)
      // + turbulence->divDevReff(U)
      - 0.5*(
              fvm::laplacian(turbulence->nuEff(), U) 
            + fvc::laplacian(turbulence->nuEff(), U.oldTime())
        )
        - fvc::div(turbulence->nuEff() * dev2((fvc::grad(U))().T()))
     ==
        - gh*fvc::grad(rhok)
        + fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve (UEqn == - fvc::grad(p_rgh) );

    fvOptions.correct(U);

