// Initialise fluid field pointer lists
PtrList<rhoReactionThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> cFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
List<scalar> rho0_value(fluidRegions.size());
List<scalar> mu_u_value(fluidRegions.size());
List<scalar> Tref_value(fluidRegions.size());
List<scalar> pref_value(fluidRegions.size());
List<scalar> ScT(fluidRegions.size());
List<scalar> X_H2_0(fluidRegions.size());
List<scalar> X_H2_99(fluidRegions.size());
List<dimensionedScalar> Y_H2_0(fluidRegions.size());
List<dimensionedScalar> Y_H2_99(fluidRegions.size());
List<scalar> X_H2O(fluidRegions.size());
List<scalar> S_L0_value(fluidRegions.size());
List<scalar> Le(fluidRegions.size());
List<Switch> ignite(fluidRegions.size());
List<Switch> visciousDis(fluidRegions.size());
PtrList<volScalarField> S_LFluid(fluidRegions.size());
PtrList<volScalarField> S_TFluid(fluidRegions.size());
PtrList<volScalarField> cSourceFluid(fluidRegions.size());
List<string> Correlation(fluidRegions.size());
List<string> CorrelationL(fluidRegions.size());
List<scalar> AC(fluidRegions.size());
List<scalar> alpha_u_value(fluidRegions.size());

#include <ctime>
std::time_t t = std::time(nullptr);
char mbstr[100];
std::strftime(mbstr, sizeof(mbstr), "%F-%H-%M-%S", std::localtime(&t));
OFstream run_info(runTime.rootPath()/runTime.globalCaseName()/"run_info."+mbstr+".txt");
run_info << "flameFoam-"+VERSION << endl;

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    run_info << "Fluid region " << fluidRegions[i].name() << " info: " << endl;
    run_info << "Mesh size: " << returnReduce(fluidRegions[i].cells().size(), sumOp<label>()) << endl;
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set(i, rhoReactionThermo::New(fluidRegions[i]).ptr());

    thermoFluid[i].validate(args.executable(), "h", "e");
    const scalar X_H2_0_c(thermoFluid[i].lookupOrDefault<scalar>("X_H2_0", 2));
    const scalar X_H2O_c(thermoFluid[i].lookupOrDefault<scalar>("X_H2O", 2));

    Info<< "    Adding to cFluid\n" << endl;
    cFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "c",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    run_info << cFluid[i].average() << endl;

    thermoFluid[i].composition().Y(0) = cFluid[i];
    volScalarField& Ydefault = thermoFluid[i].composition().Y(1);
    thermoFluid[i].correct();

    if ((cFluid[i].name() != "c") or (Ydefault.name() != "Ydefault"))
    {
        FatalErrorInFunction
        << "Wrong order of mixture components." <<  endl
        << exit(FatalError);
    }

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );

    dimensionedScalar rho0 = rhoFluid[i].average();
    rho0_value[i] = rho0.value();
    run_info << "rho0 : " << rho0 << endl;

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimLength, 0)
        )
    );

    dimensionedScalar ghRef(- mag(gFluid[i])*hRefFluid[i]);

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    run_info << "Turbulence model: " << turbulenceFluid[i].type() << endl;

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Adding to radiationFluid\n" << endl;
    radiation.set
    (
        i,
        radiationModel::New(thermoFluid[i].T())
    );

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                thermoFluid[i].p().dimensions()/dimTime,
                0
            )
        )
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulenceFluid[i].validate();

    run_info << "Average values of initial T and p:\n\t" << thermoFluid[i].T().average() << "\n\t" << thermoFluid[i].p().average() << endl;

    dimensionedScalar mu_u = average(thermoFluid[i].mu());
    mu_u_value[i] = mu_u.value();
    run_info << "mu_u : " << mu_u << endl;

    #include "initCombustion.H"
}
