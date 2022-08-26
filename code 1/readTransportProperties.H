    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Foam::Switch constMu (transportProperties.lookup("constantMu"));
    Foam::Switch constKt (transportProperties.lookup("constantKt"));

    dimensionedScalar nu0
    (
        transportProperties.lookup("nu")
    );

    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        nu0
    );       

    dimensionedScalar beta_T
    (
        transportProperties.lookup("beta_T")
    );

    dimensionedScalar T0
    (
        transportProperties.lookup("T0")
    );

    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")
    );

    // permeability!
    dimensionedScalar perm
    (
        transportProperties.lookup("perm")
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        nu0*rho0
    );

    dimensionedScalar ks (transportProperties.lookup("ks"));
    dimensionedScalar kl (transportProperties.lookup("kl"));
    dimensionedScalar Cps (transportProperties.lookup("cps"));
    dimensionedScalar Cpl (transportProperties.lookup("cpl"));
    dimensionedScalar L (transportProperties.lookup("L"));
    dimensionedScalar Tl (transportProperties.lookup("Tl"));
    dimensionedScalar Ts (transportProperties.lookup("Ts"));
    dimensionedScalar C (transportProperties.lookup("C"));
