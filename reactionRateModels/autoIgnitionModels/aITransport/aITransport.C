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

#include "aITransport.H"
#include "addToRunTimeSelectionTable.H"
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

):
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
    dataTable(),
    debug_(coeffDict_.lookupOrDefault("debug", false)),
    Sct_("Sct", dimless, 0)
{
    IOdictionary thermophysicalTransportDict
    (
        IOobject
        (
            "thermophysicalTransport",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    Sct_ = thermophysicalTransportDict.lookup<scalar>("Sct");
    appendInfo("\tAutoignition estimation method: aITransport equation");

    // Create sample data
    HashTable<scalar> table1000;
    table1000.insert(word("300"), 1.5e-3);
    table1000.insert(word("400"), 1.0e-3);
    dataTable.insert(word("1000"), table1000);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoIgnitionModels::aITransport::~aITransport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::autoIgnitionModels::aITransport::correct()
{
    if (debug_)
    {
        Info << "\t\taITransport correct:" << endl;
        Info << "\t\t\tInitial average tau: "  << average(tau_).value() << endl;
    }

    // TODO: should be done by taking thermophysicalaITransport.DEff()
    const volScalarField DL(reactionRate_.muU()/(reactionRate_.rhoU()*0.7));
    const volScalarField DT(combModel_.turbulence().nut()/Sct_);
    const volScalarField DTot(DL + DT);

    volScalarField& TU = reactionRate_.TU().ref();
    forAll (mesh_.C(), celli)
    {
        const word pRounded(Foam::name(round(p_[celli]/1000)*1000));
        const word TRounded(Foam::name(round(TU[celli])));
        ADT_[celli] = dataTable[pRounded][TRounded];
    }

    fvScalarMatrix tauEqn
    (
        fvm::ddt(rho_, tau_)
      + fvm::div(phi_, tau_)
      - fvm::laplacian(rho_*DTot, tau_)
     ==
        rho_/ADT_
    );

    // Solve equation
    tauEqn.relax();
    tauEqn.solve();

    if (debug_)
    {
        Info<< "Min/max tau: " << min(tau_).value()
            << " " << max(tau_).value() << endl;
        Info << "\t\t\taITransport correct finished" << endl;
    }
}

char const *Foam::autoIgnitionModels::aITransport::getInfo()
{
    return infoString_.c_str();
}

// ************************************************************************* //
