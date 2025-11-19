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
    const combustionModel& combModel
)
:
    combModel_(combModel),
    combustionProperties_(combModel_.coeffs()),
    mesh_(combModel_.mesh()),

    // fields
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
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    p_(mesh_.lookupObject<volScalarField>("p")),
    T_(mesh_.lookupObject<volScalarField>("T")),

    // indices and switches
    yIndex_(combModel_.thermo().specieIndex(combModel_.thermo().Y("b"))),
    Tu_(combustionProperties_.lookup<Switch>("Tu")),
    debug_(combustionProperties_.lookupOrDefault<Switch>("debug", false)), // reiktų perduot iš flameFoam
    debugFields_(combustionProperties_.lookupOrDefault<Switch>("debugFields", false)), // reiktų perduot iš flameFoam

    // model constants
    p0_(mesh_.time().value()==0 ?
        average(p_)
        : dimensionedScalar("p0", dimPressure, combustionProperties_)),

    rho0_(mesh_.time().value()==0 ?
        average(combModel_.rho())
        : dimensionedScalar("rho0", dimDensity, combustionProperties_)),

    // HEff estimation
    WH2_(dimMass/dimMoles, 2.016),
    WU_(combModel_.thermo().Wi(yIndex_)),

    X_H2_0_("X_H2_0", dimless, combustionProperties_),
    Y_H2_99_("Y_H2_99", dimless, combustionProperties_),
    Y_H2_0_(X_H2_0_*WH2_/WU_),

    H0_("H0", dimEnergy/dimMass, combustionProperties_),
    HEff_((Y_H2_0_-Y_H2_99_)*H0_),

    rhoU_(rho0_*pow(p_/p0_, 1/combModel_.thermo().gamma())),
    TU_(WU_*p_/(rhoU_*constant::physicoChemical::RR)),
    muU_(combModel_.thermo().mui(yIndex_, p_, TU_)),

    alphaU_(
        IOobject(
            "alphaU",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alpha_u", dimKinematicViscosity, combustionProperties_)),
    calculateAlphaU_(combustionProperties_.lookup<scalar>("alpha_u") <= 0.0)

{
    Info << "flameFoam reactionRate object initialized" << endl;
    Info << "Unburnt mixture temperature correction: ";
    Info << Tu_.asText() << endl;
    appendInfo(std::string("Unburnt mixture temperature correction: ") + Tu_.asText());
    appendInfo(std::string("Unburnt mixture thermal diffusivity correction: ") + calculateAlphaU_.asText());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRate::~reactionRate()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionRate::correctUnburntProperties()
{
    rhoU_ = rho0_*pow(p_/p0_, 1/combModel_.thermo().gamma());

    if (Tu_)
    {
        TU_ = WU_*p_/(rhoU_*constant::physicoChemical::RR);
    }
    else
    {
        TU_ = T_;
    }

    muU_ = combModel_.thermo().mui(yIndex_, p_, TU_);

    if (calculateAlphaU_)
    {
        alphaU_ =
            combModel_.thermo().kappai(yIndex_, p_, TU_)
            /(rhoU_*combModel_.thermo().Cpi(yIndex_, p_, TU_));
    }
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
    volScalarField& c = const_cast<volScalarField&>(combModel_.thermo().Y("c"));

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


// ************************************************************************* //
