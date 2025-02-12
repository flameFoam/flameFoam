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

#include "wFTransport.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wrinklingFactorModels
{
    defineTypeNameAndDebug(wFTransport, 0);
    addToRunTimeSelectionTable
    (
        wrinklingFactor,
        wFTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::wFTransport::wFTransport
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict

):
    wrinklingFactor(modelType, reactRate, dict),
    laminarCorrelation_
    (
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

    Le_("Le", dimless, this->coeffDict_),
    rho_(combModel_.rho()),
    phi_(mesh_.lookupObject<surfaceScalarField>("phi")),
    p_(mesh_.lookupObject<volScalarField>("p")),
    p0_("p0", dimPressure, this->coeffDict_.lookupOrDefault<scalar>("p0", 101325.0)),
    debug_(coeffDict_.lookupOrDefault("debug", false)),
    Sct_("Sct", dimless, 0)
{
    IOdictionary thermophysicalwFTransportDict
    (
        IOobject
        (
            "thermophysicalwFTransport",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    Sct_ = thermophysicalwFTransportDict.lookup<scalar>("Sct");

    appendInfo("\tWrinkling factor estimation method: wFTransport equation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wrinklingFactorModels::wFTransport::~wFTransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::wrinklingFactorModels::wFTransport::correct()
{
    if (debug_)
    {
        Info << "\t\twFTransport correct:" << endl;
        Info << "\t\t\tInitial average TBV: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();

    const volScalarField uPrime(sqrt((2.0/3.0)*combModel_.turbulence().k()));
    const volScalarField tauEta(sqrt(reactionRate_.nuU()/reactionRate_.saneEpsilon()));

    // TODO: also used in Bradley and elsewhere, probably should be in the reactionRate
    const volScalarField lT(pow(combModel_.turbulence().k(), 1.5)/reactionRate_.saneEpsilon());
    const volScalarField ReT(uPrime*lT/reactionRate_.nuU());

    // TODO: should be done by taking thermophysicalwFTransport.DEff()
    const volScalarField DL(reactionRate_.muU()/(reactionRate_.rhoU()*0.7));
    const volScalarField DT(combModel_.turbulence().nut()/Sct_);
    const volScalarField DTot(DL + DT);

    const volScalarField XiEq
    (
        scalar(1)
        + 0.46/Le_*pow(ReT, 0.25)*pow(uPrime/laminarCorrelation_->burningVelocity(), 0.3)*pow(p_/p0_, 0.2)
    );

    const volScalarField G(0.28/tauEta);
    const volScalarField R(G*XiEq/(XiEq - scalar(1)));

    // Create Xi equation
    fvScalarMatrix XiEqn
    (
        fvm::ddt(rho_, Xi_)
      + fvm::div(phi_, Xi_)  // TODO: Xi flux?
      - fvm::laplacian(rho_*DTot, Xi_)
     ==
        rho_*G*Xi_
      - rho_*R*(Xi_-scalar(1))
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
        Info << "\t\t\twFTransport correct finished" << endl;
    }
}

char const *Foam::wrinklingFactorModels::wFTransport::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}

// ************************************************************************* //
