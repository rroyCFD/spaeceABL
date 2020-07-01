// Compute the full pressure and adjust the pressure level.
// rhok is updated in TEqn
Boussinesq.updateBackgroundPressure();

p = p_rgh + Boussinesq.backgroundPressure();

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

// Shall the p, p_rgh BC be updated?