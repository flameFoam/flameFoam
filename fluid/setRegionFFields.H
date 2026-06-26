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

volScalarField& cLam = cLamFluid[i];
volScalarField& S_L = S_LFluid[i];
volScalarField& S_T = S_TFluid[i];
volScalarField& G = GFluid[i];
volScalarField& cSource = cSourceFluid[i];
volScalarField& Xi = XiFluid[i];
volScalarField& dFactor = dFactorFluid[i];
volScalarField& Pe_Q = Pe_QFluid[i];
volScalarField& Pe = PeFluid[i];
volScalarField& D_Q = D_QFluid[i];
volScalarField& quenchAlpha = quenchAlphaFluid[i];
volScalarField& I_QR = I_QRFluid[i];
volScalarField& wallDistance = wallDistanceFluid[i];

const scalar Pe_max(Pe_max_value[i]);
const scalar quenchBeta(quenchBeta_value[i]);
const dimensionedScalar rho0("rho0", dimDensity, rho0_value[i]);
const dimensionedScalar p0("p0", dimPressure, p0_value[i]);
const dimensionedScalar mu_u("mu_u", dimDynamicViscosity, mu_u_value[i]);
const dimensionedScalar W_u("W_u", dimMass/dimMoles, W_u_value[i]);
const dimensionedScalar Tref("Tref", dimTemperature, Tref_value[i]);
const dimensionedScalar pref("pref", dimPressure, pref_value[i]);
const dimensionedScalar gcr("gcr", dimless/dimTime, gcr_value[i]);
const dimensionedScalar H0("H0", dimEnergy/dimMass, 1.1996e+8);
const dimensionedScalar alpha_u("alpha_u", dimViscosity, alpha_u_value[i]);
const dimensionedScalar ACalpha("ACalpha", pow(dimViscosity, -0.25), ACalpha_value[i]);
const dimensionedScalar k_SMALL("k_SMALL", pow(dimVelocity, 2), SMALL);
const scalar LBV_param = LBV_param_value[i];
