const fvMesh& mesh = fluidRegions[i];

rhoReactionThermo& thermo = thermoFluid[i];
thermo.validate(args.executable(), "h", "e");

volScalarField& c = cFluid[i];
volScalarField& Ydefault = thermo.composition().Y(1);

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];

compressible::turbulenceModel& turbulence = turbulenceFluid[i];
volScalarField& K = KFluid[i];
volScalarField& dpdt = dpdtFluid[i];

volScalarField& p = thermo.p();
volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField& p_rgh = p_rghFluid[i];

const dimensionedVector& g = gFluid[i];
const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

radiationModel& rad = radiation[i];

IOMRFZoneList& MRF = MRFfluid[i];
fv::options& fvOptions = fluidFvOptions[i];

pimpleNoLoopControl& pimple = pimples.pimple(i);

volScalarField& S_L = S_LFluid[i];
volScalarField& S_T = S_TFluid[i];
volScalarField& cSource = cSourceFluid[i];

const dimensionedScalar rho0("rho0", dimDensity, rho0_value[i]);
const dimensionedScalar mu_u("mu_u", dimDynamicViscosity, mu_u_value[i]);
const dimensionedScalar Tref("Tref", dimTemperature, Tref_value[i]);
const dimensionedScalar pref("pref", dimPressure, pref_value[i]);
const dimensionedScalar S_L0("S_L0", dimVelocity, S_L0_value[i]);
const dimensionedScalar alpha_u("alpha_u", dimensionSet(0,-0.5,0.25,0,0,0,0), alpha_u_value[i]);
const dimensionedScalar H0("H0", dimEnergy/dimMass, 1.418e+8);
