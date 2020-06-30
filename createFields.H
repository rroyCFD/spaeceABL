// create solution control object
spaeceControl spaece(mesh, "SPAeCE");
Info << spaece.algorithmName() << endl;

Info<< "Creating/Reading thermophysical properties\n" << endl;

Info<< "Creating/Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating/Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating/Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

#include "readTransportProperties.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Create Coriolis force object [m/s^2].
Info << "Creating Coriolis force" << endl;
CoriolisForce Coriolis(U);


// Create external source term objects [m/s^2 and K/s].
Info << "Creating and reading volumetric source terms for momentum and temperature" << endl;
DrivingForce<vector> momentumSourceTerm("momentum",U);
DrivingForce<scalar> temperatureSourceTerm("temperature",T);


// // Kinematic density for buoyancy force
// volScalarField rhok
// (
//     IOobject
//     (
//         "rhok",
//         runTime.timeName(),
//         mesh
//     ),
//     1.0 - beta*(T - TRef)
// );

// // kinematic turbulent thermal thermal conductivity m2/s
// Info<< "Creating/Reading field kappat\n" << endl;
// volScalarField kappat
// (
//     IOobject
//     (
//         "kappat",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );

// Create and read the wall shear stress field [m^2/s^2].
Info << "Creating and reading the wall shear stress field, Rwall" << endl;
volSymmTensorField Rwall
(
    IOobject
    (
        "Rwall",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Create and read the wall temperature flux field [K-m/s].
Info << "Creating and reading the wall temperature flux field, qwall" << endl;
volVectorField qwall
(
    IOobject
    (
        "qwall",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// #include "readGravitationalAcceleration.H"
// #include "readhRef.H"
// #include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("", p_rgh.dimensions(), 0.0)
);

label pRefCell = 0;
scalar pRefValue = 0.0;

bool pRefOn(setRefCell (p, p_rgh, spaece.dict(), pRefCell, pRefValue ));


// Create buoyancy model that requires a reference height, 
// which is the height of the pressure reference cell.
Info << "Creating the Boussinesq buoyancy model" << endl;
vector hRef_ = vector::zero;
if (pRefCell != -1)
{
    hRef_ = mesh.C()[pRefCell];
}
reduce(hRef_,sumOp<vector>());
dimensionedVector hRef("hRef",dimLength,hRef_);
    
buoyancyModel Boussinesq(T, TRef, hRef);

// Adjust p and p_rgh with reference height and background pressure
#include "adjustPressureLevel.H"


#include "createpp.H"
mesh.setFluxRequired(pp.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
#include "createFvOptions.H"

#include "createRhieChow.H"
