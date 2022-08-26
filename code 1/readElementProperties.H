    Info<< "Reading elementProperties\n" << endl;

    IOdictionary elementProperties
    (
        IOobject
        (
            "elementProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

//      Universal gas constant (J/kg.mole.K)
    dimensionedScalar R
    (
        elementProperties.lookup("R")
    );


//      Initial concentration of oxygen (% wt)
    dimensionedScalar a_O (elementProperties.lookup("a_O_0") );

//      A_O (N/m.K) := activity => dsigmadT0
    dimensionedScalar dsigmadT0 (elementProperties.lookup("dsigmadT0"));

//      Surface excess at saturation (kg mole/m2)
    dimensionedScalar Ss_O (elementProperties.lookup("Ss_O"));

//      k1_O (dimLess)
    dimensionedScalar k1_O (elementProperties.lookup("k1_O"));

//      entropy of segregation (J/kg)
    dimensionedScalar DH_O (elementProperties.lookup("DH_O"));

//     Concentration of sulphur (% wt)
    dimensionedScalar a_S (elementProperties.lookup("a_S_0") );

//      Surface excess at saturation (kg mole/m2)
    dimensionedScalar Ss_S (elementProperties.lookup("Ss_S"));

//      k1_O (dimLess)
    dimensionedScalar k1_S (elementProperties.lookup("k1_S"));

//      entropy of segregation (J/kg)
    dimensionedScalar DH_S (elementProperties.lookup("DH_S"));

