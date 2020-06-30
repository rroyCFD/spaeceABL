{
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    // surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

    surfaceScalarField phiHbyA ("phiHbyA", fvc::flux(HbyA) );
    // + MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
    // + phig

    MRF.makeRelative(phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAU, MRF);

    // Apply the same gradient to pp following the coresponding p_rgh BCs
    forAll(p_rgh.boundaryField(), patchi)
    {
        if(isType<fixedFluxPressureFvPatchScalarField>
            (p_rgh.boundaryField()[patchi]))
        {
            fixedFluxPressureFvPatchScalarField& p_rghGradPatch =
                refCast<fixedFluxPressureFvPatchScalarField>
                    (p_rgh.boundaryFieldRef()[patchi]);

            fixedGradientFvPatchScalarField& ppGradPatch =
                refCast<fixedGradientFvPatchScalarField>
                    (pp.boundaryFieldRef()[patchi]);

            ppGradPatch.gradient() = p_rghGradPatch.gradient();

            Info << "maxGrad: " << gMax(p_rghGradPatch.gradient()) << endl;
        }
    }
}