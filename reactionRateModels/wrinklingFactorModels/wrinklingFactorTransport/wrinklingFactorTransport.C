/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2025 Lithuanian Energy Institute

 -------------------------------------------------------------------------------
License
    This file is part of flameFoam, derivative work of OpenFOAM.

    flameFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    flameFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    <http://www.gnu.org/licenses/> for more details.

Disclaimer
    flameFoam is not approved or endorsed by neither the OpenFOAM Foundation
    Limited nor OpenCFD Limited.

\*---------------------------------------------------------------------------*/

#include "wrinklingFactorTransport.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcLaplacian.H"
#include "laminarBurningVelocity.H"
#include "fvmSup.H"
#include "fvmDiv.H"
#include "fvcDdt.H"
#include "XiFluid.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wrinklingFactorModels
{
    defineTypeNameAndDebug(wrinklingFactorTransport, 0);
    addToRunTimeSelectionTable
    (
        wrinklingFactor,
        wrinklingFactorTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::wrinklingFactorTransport::wrinklingFactorTransport
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
    
):
    wrinklingFactor(modelType, reactRate, dict),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            reactRate,
            combModel_.coeffs()
        )
    ),
       Xi_
    (
        IOobject
        (
            "Xi", 
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Xi", dimless, Zero)
    ),
    SuMin_("SuMin", dimVelocity, dict.lookupOrDefault<scalar>("SuMin", 0.01)),
    SuMax_("SuMax", dimVelocity, dict.lookupOrDefault<scalar>("SuMax", 4.0)),
    XiCoef_(dict.lookupOrDefault<scalar>("XiCoef", 0.5)),
    XiShapeCoef_(dict.lookupOrDefault<scalar>("XiShapeCoef", 0.5)),
    uPrimeCoef_(dict.lookupOrDefault<scalar>("uPrimeCoef", 1.0)),
    sigmaExt_("sigmaExt", dimless/dimTime, dict.lookupOrDefault<scalar>("sigmaExt", 500)),
    Su
    (
        IOobject
        (
            "Su",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    b_(&combModel_.thermo().Y("b")),
    debug_(coeffDict_.lookupOrDefault("debug", false))
{
    appendInfo("\tWrinkling factor estimation method: wrinklingFactorTransport correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::wrinklingFactorTransport::~wrinklingFactorTransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::wrinklingFactorModels::wrinklingFactorTransport::correct()
{
     if (debug_)
    {
        Info << "\t\tCharlette correct:" << endl;
        Info << "\t\t\tInitial average TBV: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();

   //rho
   const volScalarField epsilon
    (
        pow(uPrimeCoef_, 3)*combModel_.turbulence().epsilon()
    );

    const volScalarField up(uPrimeCoef_*sqrt((2.0/3.0)*combModel_.turbulence().k()));
    const volScalarField tauEta(sqrt(reactionRate_.muU()/(reactionRate_.rhoU()*epsilon)));
     const volScalarField Reta
    (
        up
      / (
            sqrt(epsilon*tauEta)
          + dimensionedScalar(up.dimensions(), 1e-8)
        )
    );
   const volScalarField DL(reactionRate_.muU()/(reactionRate_.rhoU()*0.7));
   const volScalarField DT(combModel_.turbulence().nut()/0.7);
   const volScalarField DTot(DL + DT);
   const volScalarField XiEqStar
        (
            scalar(1.001) + XiCoef_*sqrt(up/(Su + SuMin_))*Reta
        );

    const volScalarField XiEq
    (
        scalar(1.001)
        + (
            scalar(1)
            + (2*XiShapeCoef_)
            *(scalar(0.5) - min(max(*b_, scalar(0)), scalar(1)))
        )*(XiEqStar - scalar(1.001))
    );

   const volScalarField Gstar(0.28/tauEta);
   const volScalarField R(Gstar*XiEqStar/(XiEqStar - scalar(1)));
   const volScalarField G(R*(XiEq - scalar(1.001))/XiEq);

    // Create Xi equation
    fvScalarMatrix XiEqn
    (
        fvm::ddt(reactionRate_.rhoU(), Xi_)                         
    //   + fvm::div(reactionRate_.rhoU(), Xi_)  //this is commented out because without Uj, the velocity, it's just a copy of the first line                   
      - fvm::laplacian(reactionRate_.rhoU()*DTot, Xi_)     
     ==
        reactionRate_.rhoU()*Gstar*Xi_                     
      - reactionRate_.rhoU()*R*(Xi_-scalar(1))                                   
    );

    // Solve equation
    XiEqn.relax();
    XiEqn.solve();

    // Bound Xi_ to be >= 1
    Xi_.max(1.0);
    sTurbulent_ = Xi_*laminarCorrelation_->burningVelocity();

    if (debug_)
    {
        Info<< "Min/max Xi: " << min(Xi_).value() 
            << " " << max(Xi_).value() << endl;
        Info << "\t\t\twrinklingFactorTransport correct finished" << endl;
    }
}

char const *Foam::wrinklingFactorModels::wrinklingFactorTransport::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}

// ************************************************************************* //
