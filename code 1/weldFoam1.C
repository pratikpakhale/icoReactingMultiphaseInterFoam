/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Anton Kidess
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This solver is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    The solver is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    See <http://www.gnu.org/licenses/>.

Application
    mhdWeld

Description
    Single phase laser/arc welding solver
	Author: Anton Kidess

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "Time.H"
#include "fvCFD.H"
#include "OSspecific.H"
#include "pimpleControl.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
//#   include "readGravitationalAcceleration.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
 //#include "readTimeControls.H"
#   include "initContinuityErrs.H"

    #include "postProcess.H"


    #include "createControl.H"
    #include "createTimeControls.H"



    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"


	//pimpleControl pimple(mesh);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    int maxloop = 0;

    while (runTime.run())
    {
       #include "readTimeControls.H"		
		#include "CourantNo.H"
		#include "setDeltaT.H"

        runTime++;

		Info<< "Time = " << runTime.timeName() << nl << endl;

        // temperature subcycling, need to remember lfPreloop for time derivative
        volScalarField lfPRELOOP = lf;

		//only attempt to solve velocity and pressure if the domain is not entirely solid
        while (max(T).value() >= Ts.value() && pimple.loop())
        {
            // Update solid phase dampening term
            A = -C*sqr(scalar(1)-lf)/(pow(lf,scalar(3))+ q_small)/rho0;

			// Solve the momentum equation
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              - fvm::laplacian(nu, U)
              - fvm::Sp(A, U)
            );
            UEqn.relax();

			if (pimple.momentumPredictor())
			{
				solve
				(
				    UEqn
				 ==
				    fvc::reconstruct
				    (
				        (
				          - ghf*fvc::snGrad(rhok)
				          - fvc::snGrad(p_rgh)
				        )*mesh.magSf()
				    )
				);
			}


            // --- PIMPLE loop
            while (pimple.correct())
            {
               
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U, phi)
      + phig
    );

    //MRF.makeRelative(phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

#include "continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }          
            } // pimple.correct
            // Temperature solver
            #include "TEqn.H"
        } //pimple.loop
        if (max(T) < Ts) {
            #include "TEqn.H"
        }

        runTime.write();
        Info<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
