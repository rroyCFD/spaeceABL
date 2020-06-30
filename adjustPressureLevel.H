// Compute the full pressure and adjust the pressure level.
// rhok is updated in TEqn, which solved after pressure
Boussinesq.updateBackgroundPressure();

p = p_rgh + Boussinesq.backgroundPressure();

// if (p_rgh.needReference())
if (pRefOn)
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );

    p_rgh = p - Boussinesq.backgroundPressure();
}

