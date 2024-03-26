/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "reactionRate.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reactionRate, 0);
    defineRunTimeSelectionTable(reactionRate, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRate::reactionRate
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
:
    coeffDict_(dict.optionalSubDict(modelType + "Coeffs")),
    mesh_(mesh),
    combModel_(combModel),
    mixture_(combModel_.thermo().composition()),
    cSource_
    (
        IOobject
        (
            "cSource",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cSource", dimDensity/dimTime, Zero)
    ),
    p_(mesh.lookupObject<volScalarField>("p")),
    // HEff estimation
    molarH2_(dimensionedScalar("molarH2", dimMass/dimMoles, 0.002)),
    X_H2_0_(coeffDict_.lookup<scalar>("X_H2_0")),
    Y_H2_0_(dimensionedScalar("Y_H2_0",
        molarH2_*X_H2_0_*average(p_)/(average(combModel.rho())*constant::physicoChemical::R*average(combModel.thermo().T()))
    ).value()),
    Y_H2_99_(0),
    H0_(coeffDict_.lookup<scalar>("H0")),
    HEff_(dimensionedScalar("Heff", dimEnergy/dimMass, (Y_H2_0_-Y_H2_99_)*H0_)),
    // model constants
    yIndex_(mixture_.index(mixture_.Y("b"))),
    WU_(dimensionedScalar("WU", dimMass/dimMoles, mixture_.Wi(yIndex_))),
    p0_(dimensionedScalar("p0", dimPressure, average(combModel_.thermo().p()).value())),
    rho0_(dimensionedScalar("rho0", dimDensity, average(combModel_.rho()).value())),
    // debug printout switch
    debug_(coeffDict_.lookupOrDefault("debug", false)) // reiktų perduot iš flameFoam
{
    Info << "flameFoam reactionRate object initialized" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRate::~reactionRate()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::reactionRate::TU() const
{
    return WU_*p_/(rhoU()*constant::physicoChemical::R);
}

Foam::tmp<Foam::volScalarField>
Foam::reactionRate::rhoU() const
{
    return rho0_*pow(p_/p0_, 1/combModel_.thermo().gamma());
}

Foam::tmp<Foam::volScalarField>
Foam::reactionRate::muU() const
{
    return mixture_.mu(yIndex_, p_, TU());
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::reactionRate::R(const label speciei) const
{
    if (debug_)
    {
        Info << "\treactionRate R(" << speciei << ") " << endl;
    }
    return cSource_;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::reactionRate::R(volScalarField& Y) const
{
    if (debug_)
    {
        Info << "\treactionRate R(" << Y.name() << ") " << endl;
    }

    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();
    Su += cSource_;

    return tSu;
}

Foam::tmp<Foam::volScalarField>
Foam::reactionRate::Qdot() const
{
    volScalarField& c = const_cast<volScalarField&>(mixture_.Y("c"));

    if (debug_)
    {
        Info << "\treactionRate Qdot:" << endl;
        Info << "\t\tmin/avg/max c: " << min(c).value() << " " << average(c).value() << " " << max(c).value() << endl;
    }
    tmp<volScalarField> hSource = volScalarField::New
    (
        combModel_.thermo().phasePropertyName(typedName("Qdot")),
        cSource_*HEff_*min(mag(min(c, scalar(1))-c.oldTime())/max(mag(c-c.oldTime()), VSMALL), scalar(1)) // renormalized to max c = 1
    );
    c.min(1);
    if (debug_)
    {
        volScalarField& hSourceOut = hSource.ref();
        Info << "\t\tObtained min/avg/max Qdot: " << min(hSourceOut).value() << " " << average(hSourceOut).value() << " " << max(hSourceOut).value() << endl;
        Info << "\t\tNormalized min/avg/max c: " << min(c).value() << " " << average(c).value() << " " << max(c).value() << endl;
        hSourceOut.write();
    }
    return hSource;
};



// NEAIŠKU AR REIKIA, GAL PRAVERS
// bool Foam::reactionRate::read(const dictionary& dict)
// {
//     dict.lookup("fuel") >> fuel_;
//
//     return true;
// }

// ************************************************************************* //
