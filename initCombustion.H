const dimensionedScalar molarH2 = dimensionedScalar("molarH2",dimensionSet(1,0,0,0,-1,0,0), 0.002);
const dimensionedScalar Tref("Tref", dimTemperature, 298);
const dimensionedScalar pref("pref", dimPressure, 100000);
const dimensionedScalar gcr("gcr", dimensionSet(0,0,-1,0,0,0,0), 10000.0);

Info<< "Reading combustion properties\n" << endl;
IOdictionary combustionProperties
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        fluidRegions[i],
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

X_H2_0[i] = combustionProperties.lookup<scalar>("X_H2_0");
run_info << "Initial hydrogen concentration: " << X_H2_0[i] << endl;

Y_H2_0[i] = dimensionedScalar("Y_H2_0", molarH2*X_H2_0[i]*thermoFluid[i].p().average()/(rhoFluid[i].average()*constant::physicoChemical::R*thermoFluid[i].T().average()));
Y_H2_99[i] = combustionProperties.lookup<scalar>("Y_H2_99");
run_info << "Initial hydrogen mass fraction: " << Y_H2_0[i] << endl;
run_info << "Final hydrogen mass fraction: " << Y_H2_99[i] << endl;

X_H2O[i] = combustionProperties.lookup<scalar>("X_H2O");
run_info << "Initial steam concentration: " << X_H2O[i] << endl;

if ((X_H2_0_c>1) || (X_H2O_c>1))
{
    WarningInFunction
    << "Hydrogen and/or steam concentrations not present in thermophysicalProperties." << nl << endl;
    run_info << "Hydrogen and/or steam concentrations not present in thermophysicalProperties. " << endl;

}
else
{
    if ((X_H2_0_c!=X_H2_0[i]) || (X_H2O_c!=X_H2O[i]))
    {
        FatalErrorInFunction
        << "Hyrogen and/or steam concentrations are not equal in combustionProperties and thermophysicalProperties." <<  endl
        << exit(FatalError);
    }
}

model[i] = combustionProperties.lookup<word>("model");
run_info << "Model: " << model[i] << endl;

ER[i] = 0.705*X_H2_0[i]/(0.295*(1-X_H2_0[i]-X_H2O[i]));

CorrelationL[i] = combustionProperties.lookup<word>("LBV_prediction");
if (CorrelationL[i] == "none")
{
    LBV_param_value[i] = combustionProperties.lookup<scalar>("S_L0");
    run_info << "Constant LBV used: " << LBV_param_value[i] << endl;
}
else
{
    if (CorrelationL[i] == "ANN")
    {
        Info << "ANN used, calculated equivalence ratio: " << ER[i] << endl;
        run_info << "ANN predicted LBV used, calculated equivalence ratio: " << ER[i] << endl;
        LBV_param_value[i] = ER[i];
    }
    else if (CorrelationL[i] == "Malet")
    {
        LBV_param_value[i] = 1.44*ER[i]*ER[i]+1.07*ER[i]-0.29;
        Info << "Malet's correlation used for the referential laminar velocity. Calculated equivalence ratio: " << ER[i] << ", calculated S_L0: " << LBV_param_value[i] << " m/s" << endl;
        run_info << "Malet predicted LBV used, calculated equivalence ratio: " << ER[i] << ", calculated S_L0: " << LBV_param_value[i] << " m/s" << endl;
    }
    else
    {
        FatalErrorInFunction << "Unknown LBV prediction method" << endl << exit(FatalError);
    }
}

Tref_value[i] = Tref.value();
pref_value[i] = pref.value();
gcr_value[i] = gcr.value();

Le[i] = combustionProperties.lookup<scalar>("Le");
run_info << "Le: " << Le[i] << endl;

Switch ignition(combustionProperties.lookup("ignite"));
ignite[i] = ignition;
ignX[i] = combustionProperties.lookup<scalar>("ignX");
ignY[i] = combustionProperties.lookup<scalar>("ignY");
ignZ[i] = combustionProperties.lookup<scalar>("ignZ");
ignR2[i] = pow(combustionProperties.lookup<scalar>("ignitionRadius"), 2);


S_LFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "S_L",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i],
        dimensionedScalar(dimVelocity, 0)
    )
);

S_TFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "S_T",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        S_LFluid[i]
    )
);

cSourceFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "cSource",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho0*S_TFluid[i]*mag(fvc::grad(cFluid[i]))
    )
);

turbulenceType[i] = turbulenceFluid[i].lookup<word>("simulationType");
Correlation[i] = combustionProperties.lookup<word>("correlation");

if (
    (
        turbulenceType[i] == "RAS" &&
        (
            Correlation[i] == "Zimont" ||
            Correlation[i] == "Bradley" ||
            Correlation[i] == "Bray"
        ) &&
        (
            model[i] == "TFC" ||
            model[i] == "ETFC"
        )
    )
    ||
    (
        turbulenceType[i] == "LES" &&
        (
            Correlation[i] == "PitschDuchamp" ||
            Correlation[i] == "Charlette"
        ) &&
        model[i] == "FSD"
    )
    ||
    turbulenceType[i] == "laminar"
)
{
    Info << "Turbulence model: " << turbulenceType[i] << " " << turbulenceFluid[i].type() << endl;
    run_info << "Turbulence model: " << turbulenceType[i] << " " << turbulenceFluid[i].type() << endl;
    Info << "Correlation: " << Correlation[i] << endl;
    run_info << "Correlation: " << Correlation[i] << endl;
}
else
{
    FatalErrorInFunction << "Incompatible turbulence model " << turbulenceType[i]  << " and turbulent flame speed model " << Correlation[i] << " used." << endl << exit(FatalError);
}

XiFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "Xi",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (turbulenceType[i] == "LES") ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
    )
);

quenching[i] = combustionProperties.lookup<Switch>("quenching");
if (quenching[i])
{
    run_info << "Quenching activated" << endl;
}

if (turbulenceType[i] == "RAS")
{
    volScalarField sigma = 0.26*Foam::log
    (
        0.37*Foam::pow
        (
            Foam::sqrt
            (
                turbulenceFluid[i].k()
            )
            /
            Foam::pow
            (
                thermoFluid[i].nu()*turbulenceFluid[i].epsilon(),
                0.25
            ),
            3
        )
    );
    sigma.max(SMALL);

    GFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "G",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                (quenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
            ),
            0.5*Foam::erfc
            (
                -
                (
                    Foam::log
                    (
                        15*thermoFluid[i].nu()*Foam::pow(gcr,2)
                        /
                        turbulenceFluid[i].epsilon()
                    )
                    +
                    0.5*sigma
                )
                /
                max(Foam::sqrt(2*sigma), small)
            )
        )
    );
}
else
{
    GFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "G",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            scalar(0)
        )
    );
}

TunFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "T_u",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermoFluid[i].T()
    )
);

wallQuenching[i] = combustionProperties.lookup<Switch>("wallQuenching");
if (wallQuenching[i])
{
    run_info << "Wall quenching activated" << endl;
}

wallDistanceFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "wallDistance",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (wallQuenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        wallDist(fluidRegions[i]).y()
    )
);
Pe_QFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "Pe_Q",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (wallQuenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
    )
);
PeFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "Pe",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (wallQuenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
    )
);
D_QFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "D_Q",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (wallQuenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
    )
);
quenchAlphaFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "quenchAlpha",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (wallQuenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
    )
);
I_QRFluid.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "I_QR",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            (wallQuenching[i]) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
    )
);

ScT[i] = thermophysicalTransportFluid[i].subDict(turbulenceType[i]).lookup<scalar>("Sct");
run_info << "Sc_T: " << ScT[i] << endl;

alpha_u_value[i] = combustionProperties.lookup<scalar>("alpha");
run_info << "\tUnburned alpha: " << alpha_u_value[i] << endl;
if (Correlation[i] == "Zimont")
{
    const scalar ZimontA = combustionProperties.lookup<scalar>("ZimontA");
    run_info << "\tZimont coefficient: " << ZimontA << endl;
    ACalpha_value[i] = ZimontA*Foam::pow(0.37, 0.25)*Foam::pow(alpha_u_value[i], -0.25);
}

Tu[i] = combustionProperties.lookup<Switch>("Tu");
if (Tu[i])
{
    run_info << "Fresh T estimated from pressure" << endl;
}

viscousDis[i] = combustionProperties.lookup<Switch>("viscousDis");
run_info << "Viscous dissipation: " << viscousDis[i] << endl;

Pe_max_value[i] = 40.0;
quenchBeta_value[i] = 3.7 - 2.0 * mag(ER[i]-1.0);

if (runTime.time().value() != 0.0){
    mu_u_value[i] = combustionProperties.lookup<scalar>("mu_u");
    rho0_value[i] = combustionProperties.lookup<scalar>("rho0");
    p0_value[i] = combustionProperties.lookup<scalar>("p0");
    rho0 = dimensionedScalar("rho0", dimDensity, rho0_value[i]);
    p0 = dimensionedScalar("p0", dimPressure, p0_value[i]);
    mu_u = dimensionedScalar("mu_u", dimDynamicViscosity, mu_u_value[i]);
    Y_H2_0[i] = dimensionedScalar("Y_H2_0", dimensionSet(0,0,0,0,0,0,0), combustionProperties.lookup<scalar>("Y_H2_0"));
}