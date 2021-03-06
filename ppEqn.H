/*
    rAU: inverse of the momentum equation "A" operator (the matrix diagonal) and
    rAUf: rAU interpolation from cell-centers to cell faces. 
    As a diagonal matrix, the inverse is just the diagonal matrix of the 
    reciprocal of the diagonal elements, hence the "r" in the name "rAU".
*/
volScalarField rAU(1.0/UEqn.A());
// surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

// RhieChow flux correction
if(RhieChow)
{
    tmp<surfaceVectorField> gradpDiff_f
    (
        (  (fvc::interpolate(fvc::grad(p)) & ed)
          - faceGradient.snGrad(p)
        ) * ed
    );

    phi =
        fvc::flux(U)
        + fvc::interpolate(rAU)* (gradpDiff_f() & mesh.Sf());

    gradpDiff_f.clear();
}
else
{
    phi = fvc::flux(U);
}

// ddt correction not applied as it was too diffusive
// + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi))

// MRF.makeRelative(phi);

// For cases which do not have all periodic pressure boundary.
// Adjust the balance of fluxes to obey continuity.
// normalizes flux based on mass in vs mass out balance
adjustPhi(phi, U, p);

// Update the pressure BCs to ensure flux consistency
// applicable at fixedFlux type pressure boundaries
// update fixedFluxPressure BC for p and copy gradient to pp
if(fixedFluxBC)
{
    #include "setPrghGrad.H"
}
// pressure corrector loop (with Non-orthogonal correction)
// - pp shall not set to zero and pp at old time is a better approximation than
//   0; reduced number of pressure corrections ~25% (from 14 to 11)
// pp.primitiveFieldRef() = 0.0;
// pp.correctBoundaryConditions();

while (spaece.correctNonOrthogonal())
{
    fvScalarMatrix ppEqn
    (
      fvm::laplacian(rAU, pp, "laplacian(p)") == (fvc::div(phi))
    );

    // no correction to reference cell pressure
    scalar ppRefValue = pRefValue - getRefCellValue(p, pRefCell);
    ppEqn.setReference(pRefCell, ppRefValue);

    ppEqn.solve(mesh.solver(pp.select(spaece.finalInnerIter())));

    if (spaece.finalNonOrthogonalIter())
    {
        phi -= ppEqn.flux();
    }
}

// calculate continuity error after final iteration
if(spaece.finalSPAeCEIter())
{
    #include "continuityErrs.H"
}

U -= rAU*fvc::grad(pp);
U.correctBoundaryConditions();
fvOptions.correct(U);

p_rgh += pp;
p_rgh.correctBoundaryConditions();

// set pressure correction field to zero, which reduces memory overhead
pp.primitiveFieldRef() = 0.0;

p  = p_rgh + rhok*gh;

// if (p_rgh.needReference())
if (pRefOn)
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    
    p_rgh = p - rhok*gh;
}
