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

Class
    Foam::turbulenceThermophysicalTransportModels::nonUnityLewisViscosityETFCDiffusivity

Description
    Non-unity-Lewis-ETFC-diffusivity based temperature gradient heat flux model
    for RAS or LES of turbulent flow. Allows independent specification of
    turbulent Prandtl and Schmidt numbers.  Unity laminar Lewis number is
    assumed.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

Usage
    \verbatim
    LES
    {
        model           nonUnityLewisViscosityETFCDiffusivity;
        Prt             0.85;
        Sct             0.7;
    }
    \endverbatim

SourceFiles
    nonUnityLewisViscosityETFCDiffusivity.C

\*---------------------------------------------------------------------------*/

#include "unityLewisEddyDiffusivity.H"

#ifndef nonUnityLewisViscosityETFCDiffusivity_H
#define nonUnityLewisViscosityETFCDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class nonUnityLewisViscosityETFCDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceThermophysicalTransportModel>
class nonUnityLewisViscosityETFCDiffusivity
:
    public unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
{

protected:
    // Protected data

        // Model coefficients

            //- Turbulent Schmidt number []
            dimensionedScalar Sct_;

            const volScalarField& DEffByRho_;


public:

    typedef typename TurbulenceThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        TurbulenceThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename TurbulenceThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("nonUnityLewisViscosityETFCDiffusivity");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        nonUnityLewisViscosityETFCDiffusivity
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~nonUnityLewisViscosityETFCDiffusivity()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const
        {
            return volScalarField::New
            (
                "DEff",
                this->DEffByRho_*this->momentumTransport().rho()
            );
        }

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff() const
        {
            return volScalarField::New
            (
                "DEff",
                this->DEffByRho_*this->momentumTransport().rho()
            );
        }

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const
        {
            return this->DEffByRho_.boundaryField()[patchi]*this->momentumTransport().rho().boundaryField()[patchi];
        }

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonUnityLewisViscosityETFCDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
