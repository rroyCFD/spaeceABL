// Calculate divergence of velocity flux and display
{
    // Compute divergence cell by cell and report statistics.  
    volScalarField divPhiMag  = mag(fvc::div(phi));

    scalar minLocalPhiContErr = min(divPhiMag).value();
    reduce(minLocalPhiContErr, minOp<scalar>());

    scalar maxLocalPhiContErr = max(divPhiMag).value();
    reduce(maxLocalPhiContErr, maxOp<scalar>());
    
    scalar avgLocalPhiContErr = divPhiMag.weightedAverage(mesh.V()).value();
    Info<< "Continuity Report:" << nl
        << "   -Local Cell Continuity Error:" << nl
        << "       minimum: " << minLocalPhiContErr << nl
        << "       maximum: " << maxLocalPhiContErr << nl
        << "       weighted mean: " << avgLocalPhiContErr << endl;
}

{
    // Compute global divergence over domain boundaries and report statistics.
    Info << "   -Boundary Flux:" << endl;

    scalar globalSumPhiBoundary(0.0);
    // scalar globalSumAreaBoundary(0.0);
    forAll(phi.boundaryField(), patchi)
    {
        // coupled BCs are always flux balanced
        if(mesh.boundary()[patchi].coupled())
        {
            continue;
        }
        // assuming non-permeable walls
        else if(isType<wallFvPatch>(mesh.boundary()[patchi]))
        {
            continue;
        }    
        else    
        {
            word boundaryName = mesh.boundary()[patchi].name();
            scalar sumPhiBoundary(0.0);
            // scalar sumAreaBoundary(0.0);
            const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
            forAll(phip,i)
            {
                sumPhiBoundary += phip[i];
                // sumAreaBoundary += mesh.boundary()[patchi].magSf()[i];
            }

            reduce(sumPhiBoundary, sumOp<scalar>());
            // reduce(sumAreaBoundary, sumOp<scalar>());

            globalSumPhiBoundary += sumPhiBoundary;
            // globalSumAreaBoundary += sumAreaBoundary;

            // Info << "       " << boundaryName << " - flux: " << sumPhiBoundary << tab << "/ area: " << sumAreaBoundary << endl;
            Info << "       " << boundaryName << " - flux: " << sumPhiBoundary << endl;
        }
    }

    // Info << "       total - flux: " << globalSumPhiBoundary << tab << "/ area: " << globalSumAreaBoundary << endl;
    Info << "       total - flux: " << globalSumPhiBoundary << endl;
}
