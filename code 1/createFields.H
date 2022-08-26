	// momentum source term 'small' constant
    dimensionedScalar q_small(transportProperties.lookup("qsmall"));

    Info<< "Reading field T\n" << endl;
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

    //store temperature gradient for Marangoni boundary condition
    volVectorField Tg
    (
        IOobject
        (
            "Tg",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(T)
    );


	Info<< "Reading field p_rgh\n" << endl;
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

    Info<< "Reading field U\n" << endl;
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

#   include "createPhi.H"

//element properties -> Oxygen, dsigmadT
#   include "readElementProperties.H"

	//liquid fraction
    volScalarField lf
    (
        IOobject
        (
            "lf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta_T*(T - T0)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


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
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    //pimple.dict(),
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

	// thermal conductivity
    volScalarField kT
    (
        IOobject
        (
            "kT",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        ks,
        lf.boundaryField().types() //zeroGradient
    );

	// surface tension variation with temperature / marangoni coefficient
    volScalarField dsigmadT 
    (
        IOobject
        (
            "dsigmadT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        -dsigmadT0,
        lf.boundaryField().types() //lf has zerogradient on all boundaries, which is suitable for dsigmadT as well
    );

	// Momentum sink coefficient
    volScalarField A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -C*sqr(1.0-lf)/(pow(lf,3.0)+q_small)/rho0
    );

Info<< "Creating field J\n" << endl;
volVectorField J
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), Foam::vector(0,0,0))
);

// latent heat source terms (implicit)
    volScalarField Sp
    (
        IOobject
        ( 
            "Sp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Sp", dimensionSet(0,2,-3,-1,0,0,0), 0)
    );

// latent heat source terms (explicit)
    volScalarField Sc
    (
        IOobject
        ( 
            "Sc",
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Sc", dimensionSet(0,2,-3,0,0,0,0), 0)
    );



