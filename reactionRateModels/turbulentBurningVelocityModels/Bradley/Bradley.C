/*---------------------------------------------------------------------------*\

 flameFoam
 Copyright (C) 2021-2024 Lithuanian Energy Institute

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

#include "Bradley.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentBurningVelocityModels
{
    defineTypeNameAndDebug(Bradley, 0);
    addToRunTimeSelectionTable
    (
        turbulentBurningVelocity,
        Bradley,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Bradley::Bradley
(
    const word modelType,
    const reactionRate& reactRate,
    const dictionary& dict
):
    turbulentBurningVelocity(modelType, reactRate, dict),
    Le_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("Le"))
{
    appendInfo("\tTBV estimation method: Bradley correlation");
    appendInfo("\t\tLe: " + name(Le_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentBurningVelocityModels::Bradley::~Bradley()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::turbulentBurningVelocityModels::Bradley::correct()
{
    if (debug_)
    {
        Info << "\t\tBradley correct:" << endl;
        Info << "\t\t\tInitial average S_T: "  << average(sTurbulent_).value() << endl;
    }

    laminarCorrelation_->correct();
    sTurbulent_ =
        1.37179015019233* //0.88*0.157^(-0.3)*(2/3)^(0.275)
        pow(laminarCorrelation_->burningVelocity(), 0.6)*pow(saneEpsilon()*reactionRate_.muU()/reactionRate_.rhoU(), -0.15)*Foam::pow(Le_, -0.3)*pow(combModel_.turbulence().k(), 0.5);

    if (debug_)
    {
        Info << "\t\t\tObtained average S_T: "  << average(sTurbulent_).value() << endl;
        Info << "\t\t\tBradley correct finished" << endl;
    }
    // dimensionedScalar omega0
    // (
    //     "omega0",
    //     dimensionSet(1, -2, -1, 0, 0, 0, 0),
    //     correlation_.omega0()
    // );
    //
    // dimensionedScalar sigmaExt
    // (
    //     "sigmaExt",
    //     dimensionSet(0, 0, -1, 0, 0, 0, 0),
    //     correlation_.sigmaExt()
    // );
    //
    // dimensionedScalar omegaMin
    // (
    //     "omegaMin",
    //     omega0.dimensions(),
    //     1e-4
    // );
    //
    // dimensionedScalar kMin
    // (
    //     "kMin",
    //     sqr(dimVelocity),
    //     small
    // );
    //
    // const compressibleMomentumTransportModel& turbulence =
    //     combModel_.turbulence();
    //
    // // Total strain
    // const volScalarField sigmaTotal
    // (
    //     sigma + alpha_*turbulence.epsilon()/(turbulence.k() + kMin)
    // );
    //
    // const volScalarField omegaInf(correlation_.omega0Sigma(sigmaTotal));
    //
    // dimensionedScalar sigma0("sigma0", sigma.dimensions(), 0.0);
    //
    // const volScalarField tau(C_*mag(sigmaTotal));
    //
    // volScalarField Rc
    // (
    //     (tau*omegaInf*(omega0 - omegaInf) + sqr(omegaMin)*sigmaExt)
    //    /(sqr(omega0 - omegaInf) + sqr(omegaMin))
    // );
    //
    // const volScalarField& rho = combModel_.rho();
    // const tmp<surfaceScalarField> tphi = combModel_.phi();
    // const surfaceScalarField& phi = tphi();
    //
    // solve
    // (
    //      fvm::ddt(rho, omega_)
    //    + fvm::div(phi, omega_)
    //   ==
    //      rho*Rc*omega0
    //    - fvm::SuSp(rho*(tau + Rc), omega_)
    // );
    //
    // omega_.min(omega0);
    // omega_.max(0.0);


}

char const *Foam::turbulentBurningVelocityModels::Bradley::getInfo()
{
    infoString_.append(laminarCorrelation_().getInfo());
    laminarCorrelation_().clearInfo();
    return infoString_.c_str();
}


// bool  Foam::turbulentBurningVelocityModels::Bradley::read
// (
//     const dictionary& dict
// )
// {
//     if (reactionRateFlameArea::read(dict))
//     {
//         coeffDict_ = dict.optionalSubDict(typeName + "Coeffs");
//         coeffDict_.lookup("C") >> C_;
//         coeffDict_.lookup("alpha") >> alpha_;
//         correlation_.read
//         (
//             coeffDict_.subDict(fuel_)
//         );
//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }

// ************************************************************************* //
