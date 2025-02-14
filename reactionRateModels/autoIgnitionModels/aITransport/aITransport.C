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
#include <dirent.h>
#include <sstream>
#include "IFstream.H"


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

      // Load the ADT data
    loadADTData();
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

void Foam::autoIgnitionModels::aITransport::loadADTData()
{
    const char* dataPath = "constant/ADT";
    DIR* dir = opendir(dataPath);
    
    if (!dir)
    {
        FatalErrorInFunction
            << "Cannot open directory " << dataPath
            << exit(FatalError);
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr)
    {
        string filename(entry->d_name);
        
        // Skip . and .. directories and non-.ADT files
        if (filename == "." || filename == ".." || 
            filename.substr(filename.length() - 4) != ".ADT")
        {
            continue;
        }

        // Use filename without .ADT extension as key
        word pKey = filename.substr(0, filename.length() - 4);
        
        // Create inner table for this pressure
        HashTable<scalar> innerTable;
        
        // Read and parse file
        string fullPath = string(dataPath) + "/" + filename;
        IFstream dataFile(fullPath);
        string line;
        
        // Skip header (first 4 lines)
        for(int j = 0; j < 4; j++)
        {
            dataFile.getLine(line);
        }
        
        // Read data lines
        while (dataFile.good())
        {
            dataFile.getLine(line);
            if (line.empty()) continue;
            
            // Parse temperature and ignition time
            scalar temp, igTime;
            std::istringstream iss(line);
            iss >> temp >> igTime;
            
            if (iss)
            {
                // Round temperature to nearest integer for table key
                word tempKey = Foam::name(round(temp));
                innerTable.insert(tempKey, igTime);
            }
        }
        
        // Add to main table
        dataTable.insert(pKey, innerTable);
    }
    
    closedir(dir);

    if (debug_)
    {
        Info<< "\nLoaded ADT data for " << dataTable.size() << " pressure values" << endl;
        
        // Print random 10 entries of main table
        Info<< "\nRandom 10 pressure files in main table:" << endl;
        label count = 0;
        forAllConstIter(HashTable<HashTable<scalar>>, dataTable, iter)
        {
            if (count++ >= 10) break;
            Info<< "  " << iter.key() << endl;
        }
        
        // Print random 10 entries of a random inner table
        if (!dataTable.empty())
        {
            const HashTable<scalar>& firstInnerTable = dataTable.begin()();
            Info<< "\nRandom 10 temperature entries for " << dataTable.begin().key() << ":" << endl;
            count = 0;
            forAllConstIter(HashTable<scalar>, firstInnerTable, innerIter)
            {
                if (count++ >= 10) break;
                Info<< "  T = " << innerIter.key() << " K, tau = " << innerIter() << " s" << endl;
            }
        }
        
        Info<< endl;
    }
}

char const *Foam::autoIgnitionModels::aITransport::getInfo()
{
    return infoString_.c_str();
}

// ************************************************************************* //
