// variables for convergence monitoring
int loopcount = 0;
scalar conv = 1;

// read subcycle settings from PISO-subdict
label LfMaxIter = 10; //(readLabel(pimple.lookup("LfMaxIter")));
scalar LfConverge (readScalar (mesh.solutionDict().subDict("PISO").lookup("LfConverge")));

while (conv > LfConverge) {

    loopcount++;
    maxloop = max(loopcount, maxloop); // has to be declared outside time loop to keep value

    volScalarField lfold = lf; // previous subcycle value

    // source term variables
    dimensionedScalar dFdT = 1 / (Tl - Ts); // depends on nature of update equation, not generally constant
    volScalarField invF = lfold * (Tl - Ts) + Ts; // invF returns liquid fraction for a given Temperature

    //update thermal transport properties if kl != ks
    if (!constKt)
    {
        //Mills, Keane
        Info << "Update kT" << endl;
        dimensionedScalar keff298 = ks;
        dimensionedScalar ks_switch ("ks_switch", dimensionSet(0,0,0,1,0,0,0), 1060.0);
        dimensionedScalar T298("T298", dimensionSet(0,0,0,1,0,0,0), 298.0);
        dimensionedScalar T775("T775", dimensionSet(0,0,0,1,0,0,0), 775.0);
        dimensionedScalar T1073("T1073", dimensionSet(0,0,0,1,0,0,0), 1073.0);
        dimensionedScalar kbase ("kbase", dimensionSet(1,1,-3,-1,0,0,0), 25.0);
        dimensionedScalar kfactor ("kfactor", dimensionSet(1,1,-3,-2,0,0,0), 0.013);
        
        volScalarField ks_T = (keff298 + (T-T298)/T775 * (kbase - keff298)) * pos(ks_switch - T) + (kbase + kfactor * (T - T1073)) * pos(T - ks_switch);
        // linear interpolation for kT
        kT = ks_T * (1.0 - lf) + kl * lf;
    } else {
        kT = ks * (1.0 - lf) + kl * lf;
    }
    
    kT.correctBoundaryConditions();
    T.correctBoundaryConditions();

    volScalarField Cp = lf * Cpl + (1-lf) * Cps;
    surfaceScalarField phiCp = fvc::interpolate(Cp) * phi;
    volScalarField dH = L * pos(Tl - T) * pos(T - Ts);

    // two source term contributions, limited to mushy zone temperatures
    Sp = - dH/runTime.deltaT() * dFdT;
    Sc = dH/runTime.deltaT() * (lfPRELOOP - lf + dFdT * invF);

    fvScalarMatrix TEqn
    (
        fvm::ddt(Cp, T)
      + fvm::div(phiCp, T)
      - fvm::laplacian(kT/rho0,T)
      - fvm::Sp(Sp, T) // latent heat 1: implicit
     ==
        Sc // latent heat 2: explicit
    );

    solve ( TEqn );

    // liquid fraction correction, convergence criterion
    conv = -1;
    lf = (lfold + (dFdT*(T-invF)));
    lf = max(lf,scalar(0));;
    lf = min(scalar(1),lf);
    conv = max(mag((lf-lfold)/(lfold+scalar(1e-3)))).value();

    if (loopcount > LfMaxIter) break;
    Info << "\tloopcount: " << loopcount << ", \tRESIDUAL: " << conv << endl;
}
Info << "\t" << loopcount << " temperature subcycles (max: " << maxloop << "). Final residual: " << conv << endl;
Info << "Min / Max T: " << min(T) << " / " << max(T) << endl;

//update auxillary fields
Tg == fvc::grad(T);
rhok = 1.0 - beta_T*(T-T0);

//Update surface tension gradient
volScalarField Kseg_S = k1_S*Foam::exp(-DH_S/(R*T));
volScalarField Kseg_O = k1_O*Foam::exp(-DH_O/(R*T));
dsigmadT = 
    (
        - dsigmadT0 
        - R*Ss_S*Foam::log(1+Kseg_S*a_S) - Kseg_S*a_S*Ss_S*DH_S/((1+Kseg_S*a_S)*T)
        - R*Ss_O*Foam::log(1+Kseg_O*a_O) - Kseg_O*a_O*Ss_O*DH_O/((1+Kseg_O*a_O)*T)
    ) * pos(T-Tl);
dsigmadT.correctBoundaryConditions();

//Update viscosity
if (!constMu)
{
    dimensionedScalar unitmu("unitmu", dimMass/(dimLength*dimTime), 1.0);
    dimensionedScalar unitT("unitT", dimTemperature, 1.0);
    nu = unitmu*exp(-2.396+(7950.0*unitT/(T+unitT)))/(1000.0*rho0) * pos(T - Tl) + nu0 * pos(Tl - T);
    mu = nu * rho0;
    Info << "Min / Max nu: " << min(nu) << " / " << max(nu) << endl;
} else {
    mu = nu * rho0;
}
U.correctBoundaryConditions();