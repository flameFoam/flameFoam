/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2026 Lithuanian Energy Institute

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

#include "aITransport.H"
#include "addToRunTimeSelectionTable.H"
#include "lookupSct.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace autoIgnitionModels
{
    defineTypeNameAndDebug(aITransport, 0);
    addToRunTimeSelectionTable
    (
        autoIgnition,
        aITransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoIgnitionModels::aITransport::aITransport
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
)
:
    autoIgnition(modelType, reactRate, dict),
    rho_(combModel_.rho()),
    phi_(mesh_.lookupObject<surfaceScalarField>("phi")),
    p_(mesh_.lookupObject<volScalarField>("p")),
    ADT_
    (
        IOobject
        (
            "ADT",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("ADT", dimTime, scalar(1))
    ),
    adtTable_
    (
        dict.lookupOrDefault<fileName>("ADTDir", fileName("constant/ADT")),
        dict.lookupOrDefault<label>("nHeaderLines", 4)
    ),
    Sct_("Sct", dimless, -1)
{
    reactionRate_.appendInfo
    (
        "\tAutoignition estimation method: aITransport equation"
    );
    reactionRate_.appendInfo
    (
        "\t\tADT table pressure samples: " + Foam::name(adtTable_.nX())
    );
    reactionRate_.appendInfo
    (
        "\t\tADT header lines skipped: "
      + Foam::name(dict.lookupOrDefault<label>("nHeaderLines", 4))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoIgnitionModels::aITransport::~aITransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::autoIgnitionModels::aITransport::correct()
{
    if (Sct_.value() < 0)
    {
        Sct_ = lookupSctFromRegistry(mesh_);
    }

    if (debug_)
    {
        Info<< "\t\taITransport correct:" << endl;
        Info<< "\t\t\tInitial average tau: " << average(tau_).value() << endl;
    }

    const volScalarField DL(reactionRate_.muU()/(reactionRate_.rhoU()*0.7));
    const volScalarField DT(combModel_.turbulence().nut()/Sct_);
    const volScalarField DTot(DL + DT);

    const volScalarField& TU = reactionRate_.TU();
    forAll(mesh_.C(), celli)
    {
        ADT_[celli] = adtTable_.interpolate(p_[celli], TU[celli]);
    }
    ADT_.correctBoundaryConditions();

    fvScalarMatrix tauEqn
    (
        fvm::ddt(rho_, tau_)
      + fvm::div(phi_, tau_)
      - fvm::laplacian(rho_*DTot, tau_)
     ==
        rho_/ADT_
    );

    tauEqn.relax();
    tauEqn.solve();

    if (debug_)
    {
        Info<< "Min/max tau: " << min(tau_).value()
            << " " << max(tau_).value() << endl;
        Info<< "\t\t\taITransport correct finished" << endl;
    }
}


// ************************************************************************* //
