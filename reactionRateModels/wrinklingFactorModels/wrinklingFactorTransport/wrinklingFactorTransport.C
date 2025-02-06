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
    U_(mesh_.lookupObject<volVectorField>("U")),
    delta_(mesh_.objectRegistry::lookupObject<volScalarField>("delta")),
    laminarCorrelation_(
        laminarBurningVelocity::New
        (
            reactRate,
            combModel_.coeffs()
        )
    ),
    c2_(0.5),
    Ck_(1.5),
    Ck_mult1_(0.245454545454545), // 27/110
    Ck_mult2_(0.648567745274438), // 4*sqrt(27/110)*18/55
    n43_(4.0/3.0),
    pi43_(Foam::pow(Foam::constant::mathematical::pi, n43_)),
    beta_(0.5)
{
    appendInfo("\tWrinkling factor estimation method: wrinklingFactorTransport correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::wrinklingFactorTransport::~wrinklingFactorTransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::wrinklingFactorModels::wrinklingFactorTransport::correct()
{
   volScalarField& Xi(Xi_);
   //rho
   //muj
   //tauEta
   //XiEqStar
   //XiEq
   //DTot = DL + DT
   //DL
   //DT
   const volScalarField Gstar(0.28/tauEta);
   const volScalarField R(Gstar*XiEqStar/(XiEqStar - scalar(1)));
   const volScalarField G(R*(XiEq - scalar(1.001))/XiEq);

    // Create Xi equation
    fvScalarMatrix XiEqn
    (
        fvm::ddt(rho, Xi)                         
      + fvm::div(phi*muj, Xi)                     
      - fvm::laplacian(rho*DTot, Xi)     
     ==
        rho*Gstar*Xi                     
      - rho*R*(Xi-scalar(1))                                   
    );

    // Solve equation
    XiEqn.relax();
    XiEqn.solve();

    // Bound Xi to be >= 1
    Xi_.max(1.0);

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
