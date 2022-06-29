const fvMesh& mesh = fluidRegions[i];

fluidReactionThermo& thermo = thermoFluid[i];
thermo.validate(args.executable(), "h", "e");

volScalarField& c = cFluid[i];
volScalarField& Ydefault = thermo.composition().Y(1);

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];

compressible::momentumTransportModel& turbulence = turbulenceFluid[i];
fluidReactionThermophysicalTransportModel& thermophysicalTransport =
        thermophysicalTransportFluid[i];

volScalarField& K = KFluid[i];
volScalarField& dpdt = dpdtFluid[i];

volScalarField& p = thermo.p();
volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField& p_rgh = p_rghFluid[i];

const dimensionedVector& g = gFluid[i];
const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
    fieldsFluid[i];

IOMRFZoneList& MRF = MRFfluid[i];
Foam::fvModels& fvModels = fvModelsFluid[i];
Foam::fvConstraints& fvConstraints = fvConstraintsFluid[i];

#include "checkRadiationModel.H"

pimpleNoLoopControl& pimple = pimples.pimple(i);

// Placeholder, will be used for moving mesh support
autoPtr<surfaceVectorField> rhoUf(nullptr);

volScalarField& S_L = S_LFluid[i];
volScalarField& S_T = S_TFluid[i];
volScalarField& cSource = cSourceFluid[i];
volScalarField& Xi = XiFluid[i];

const dimensionedScalar rho0("rho0", dimDensity, rho0_value[i]);
const dimensionedScalar mu_u("mu_u", dimDynamicViscosity, mu_u_value[i]);
const dimensionedScalar Tref("Tref", dimTemperature, Tref_value[i]);
const dimensionedScalar pref("pref", dimPressure, pref_value[i]);
const dimensionedScalar H0("H0", dimEnergy/dimMass, 1.418e+8);
const dimensionedScalar alpha_u("alpha_u", dimViscosity, alpha_u_value[i]);
const dimensionedScalar ACalpha("ACalpha", pow(dimViscosity,-0.25), ACalpha_value[i]);
const scalar LBV_param = LBV_param_value[i];
