const dimensionedScalar molarH2 = dimensionedScalar("molarH2",dimensionSet(1,0,0,0,-1,0,0), 0.002);
const dimensionedScalar Tref("Tref", dimTemperature, 298);
const dimensionedScalar pref("pref", dimPressure, 100000);

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
Y_H2_99[i] = 0;
run_info << "Initial hydrogen mass fraction: " << Y_H2_0[i] << endl;

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

CorrelationL[i] = combustionProperties.lookup<word>("LBV_prediction");
if (CorrelationL[i] == "none")
{
    LBV_param_value[i] = combustionProperties.lookup<scalar>("S_L0");
    run_info << "Constant LBV used: " << LBV_param_value[i] << endl;
}
else
{
    const scalar ER = 0.705*X_H2_0[i]/(0.295*(1-X_H2_0[i]-X_H2O[i]));
    if (CorrelationL[i] == "ANN")
    {
        Info << "ANN used, calculated equivalence ratio: " << ER << endl;
        run_info << "ANN predicted LBV used, calculated equivalence ratio: " << ER << endl;
        LBV_param_value[i] = ER;
    }
    else if (CorrelationL[i] == "Malet")
    {
        LBV_param_value[i] = 1.44*ER*ER+1.07*ER-0.29;
        Info << "Malet's correlation used for the referential laminar velocity. Calculated equivalence ratio: " << ER << ", calculated S_L0: " << LBV_param_value[i] << " m/s" << endl;
        run_info << "Malet predicted LBV used, calculated equivalence ratio: " << ER << ", calculated S_L0: " << LBV_param_value[i] << " m/s" << endl;
    }
    else
    {
        FatalErrorInFunction << "Unknown LBV prediction method" << endl << exit(FatalError);
    }
}

Tref_value[i] = Tref.value();
pref_value[i] = pref.value();

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
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i],
        scalar(0)
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

if (((turbulenceType[i] == "RAS") && ((Correlation[i] == "Zimont") || (Correlation[i] == "Bradley") || (Correlation[i] == "Bray"))) || ((turbulenceType[i] == "LES") && ((Correlation[i] == "PitschDuchamp") || (Correlation[i] == "Charlette") )) || (turbulenceType[i] == "laminar"))
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

Switch visc_switch(combustionProperties.lookup("visciousDis"));
visciousDis[i] = visc_switch;
run_info << "Viscious dissipation: " << visciousDis[i] << endl;
