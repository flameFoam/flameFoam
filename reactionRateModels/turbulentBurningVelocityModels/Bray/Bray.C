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

#include "Bray.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{
    defineTypeNameAndDebug(Bray, 0);
    addToRunTimeSelectionTable
    (
        turbulentBurningVelocity,
        Bray,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Bray::Bray
(
    const dictionary& dict,
    const reactionRate& reactRate
):
    turbulentBurningVelocity(reactRate),
    c1_(0.875 * pow(0.157, -0.392) * pow(2.0/3.0, 0.206))
{
    appendInfo("\tTBV estimation method: Bray correlation");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Bray::~Bray()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::turbulentBurningVelocityModels::Bray::correct()
{
    if (debug_)
    {
        Info << "\t\tBray correct:" << endl;
        Info << "\t\t\tInitial average S_T: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();
    sTurbulent_ =
        c1_
        *pow(combModel_.turbulence().k(), 0.5)
        *pow(laminarCorrelation_->burningVelocity(), 0.784)
        *pow(saneEpsilon()*reactionRate_.muU()/reactionRate_.rhoU(), -0.196);

    if (debug_)
    {
        Info << "\t\t\tObtained average S_T: "  << average(sTurbulent_).value() << endl;
        Info << "\t\t\tBray correct finished" << endl;
    }
}

char const *Foam::turbulentBurningVelocityModels::Bray::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}
