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

#include "ANN.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarBurningVelocityModels
{
    defineTypeNameAndDebug(ANN, 0);
    addToRunTimeSelectionTable
    (
        laminarBurningVelocity,
        ANN,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::ANN::ANN
(
    const word modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel,
    const reactionRate& reactRate
):
    laminarBurningVelocity(modelType, dict, mesh, combModel, reactRate),
    X_H2_0_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2_0")),
    X_H2O_(dict.optionalSubDict(modelType + "Coeffs").lookup<scalar>("X_H2O")),
    ER_(0.705*X_H2_0_/(0.295*(1-X_H2_0_-X_H2O_))),
    p_(mesh.lookupObject<volScalarField>("p")),

    par_(PtrList<volScalarField>(3)),

    // Weight matrices initialisation   
    W0_(scalarListList({{-3.243977308273315430e+00,
    -8.332144469022750854e-03,
    -3.480189561843872070e+00,
    -3.300335407257080078e+00,
    -1.697143077850341797e+00,
    1.478089523315429688e+01,
    -4.349520683288574219e+00}, 
    {3.029801368713378906e+00,
    1.122648711316287518e-03,
    -1.612274795770645142e-01,
    1.099239826202392578e+00,
    5.997737884521484375e+00,
    6.218534111976623535e-01,
    2.587657421827316284e-03}, 
    {-8.222029209136962891e-01,
    -1.702534849755465984e-03,
    1.868390440940856934e+00,
    -1.437355875968933105e+00,
    7.229740265756845474e-03,
    -5.841984152793884277e-01,
    -1.154016256332397461e+00}
    })),

    B0_(scalarList({
    -4.164304435253143311e-01,
    7.309035572689026594e-05,
    -5.290837287902832031e-01,
    4.410822391510009766e-01,
    2.554211765527725220e-02,
    -5.750028043985366821e-02,
    6.104745268821716309e-01
    })),

    W1_(scalarListList({
    {5.511226877570152283e-02,
    4.298480344004929066e-04,
    -8.532372117042541504e-02,
    7.387475371360778809e-01,
    1.186388079077005386e-03,
    -1.232936605811119080e-01,
    -4.893136210739612579e-03,
    -4.134190678596496582e-01,
    1.398074673488736153e-04,
    -3.192920386791229248e-01},
    {2.666092186700552702e-04, 
    2.076551900245249271e-04,
    -3.022944147232919931e-04,
    -1.304411562159657478e-04,
    2.771862782537937164e-04,
    1.071571314241737127e-04,
    4.969434812664985657e-04,
    7.113337051123380661e-05,
    2.117981784977018833e-05,
    -3.026006161235272884e-04}, 
    {1.559144351631402969e-03,
    3.275305207353085279e-04,
    -2.183801494538784027e-02,
    -4.137277901172637939e-01,
    5.053054094314575195e-01,
    -6.569375749677419662e-03,
    -3.687434364110231400e-03,
    -7.323424816131591797e-01,
    -8.242395124398171902e-04,
    -9.765876084566116333e-02},
    {6.598435044288635254e-01,
    -2.240608446300029755e-04,
    -7.942664087750017643e-04,
    1.797222648747265339e-04,
    -4.500378966331481934e-01,
    1.138506224378943443e-03,
    1.627122052013874054e-03,
    1.399135589599609375e-01,
    -6.955675780773162842e-05,
    -1.075079366564750671e-01},
    {-5.974398180842399597e-02,
    1.723054156173020601e-04,
    1.273235797882080078e+00,
    -1.731790543999522924e-04,
    -8.303441107273101807e-02,
    -3.123809695243835449e-01,
    1.995624154806137085e-01,
    1.835895180702209473e-01,
    8.467041850090026855e-01,
    3.322252035140991211e-01},
    {-4.169672320131212473e-04,
    -4.472784348763525486e-05,
    -1.608545780181884766e-01,
    2.697963640093803406e-02,
    8.455475326627492905e-04,
    -2.726765871047973633e-01,
    -1.094287753105163574e+00,
    1.316546440124511719e+00,
    -3.417208790779113770e-04,
    9.424296021461486816e-01}, 
    {4.223680496215820312e-01,
    -8.689393871463835239e-05,
    -8.681480772793292999e-03,
    -1.382667571306228638e-01,
    -8.501461893320083618e-02,
    6.597061157226562500e-01,
    6.922888159751892090e-01,
    -7.694861851632595062e-04,
    6.948460941202938557e-04,
    2.327049849554896355e-03}
    })),

    B1_(scalarList({
    4.296131730079650879e-01,
    -1.001953380182385445e-03,
    -2.080285549163818359e-01,
    2.410185933113098145e-01,
    3.027851581573486328e-01,
    5.177603363990783691e-01,
    9.652259945869445801e-02,
    5.754043534398078918e-02,
    4.636074975132942200e-02,
    -1.621754169464111328e-01
    })),

    W2_(scalarListList({
    {-4.783591139130294323e-04,
    -2.981348931789398193e-01,
    4.586713612079620361e-01,
    -1.021348163485527039e-01,
    -4.453010915312916040e-04,
    8.886648574844002724e-04,
    -5.310251116752624512e-01},
    {-1.485306711401790380e-04,
    2.168068313039839268e-04,
    -3.013435052707791328e-04,
    2.161701850127428770e-04,
    -3.544702776707708836e-04,
    4.538447537925094366e-04,
    -1.486202527303248644e-04},
    {5.931606888771057129e-01,
    4.774236381053924561e-01,
    -8.324894309043884277e-01,
    5.146421864628791809e-02,
    3.168596886098384857e-04,
    2.196232602000236511e-02,
    2.066201996058225632e-04},
    {1.064987704157829285e-01,
    -6.327702403068542480e-01,
    -4.581197351217269897e-02,
    1.660551351960748434e-04,
    3.450178191997110844e-04,
    -1.972953826189041138e-01,
    3.245791792869567871e-01},
    {-2.085151374340057373e-01,
    6.744181737303733826e-04,
    4.834557175636291504e-01,
    8.392338454723358154e-02,
    2.140337601304054260e-04,
    -5.556737184524536133e-01,
    -5.480980616994202137e-04},
    {3.176291286945343018e-01,
    -6.121918559074401855e-01,
    -6.714318878948688507e-04,
    -1.777123659849166870e-03,
    -1.299364957958459854e-04,
    4.049893934279680252e-03,
    -3.300224198028445244e-04},
    {9.933376312255859375e-01,
    1.057758199749514461e-04,
    -4.851880366913974285e-04,
    -1.376400738954544067e-01,
    -2.388842403888702393e-07,
    3.151315904688090086e-04,
    4.596516489982604980e-01},
    {3.815629184246063232e-01,
    -4.374919533729553223e-01,
    1.985776424407958984e-02,
    -1.122563481330871582e+00,
    -4.319254076108336449e-04,
    -1.456504454836249352e-03,
    -6.638678312301635742e-01},
    {2.396315187215805054e-01,
    3.806896209716796875e-01,
    3.197264741174876690e-04,
    1.720842570066452026e-01,
    3.279441734775900841e-04,
    6.297810468822717667e-04,
    3.376027047634124756e-01},
    {-1.050056442618370056e-01,
    5.293983817100524902e-01,
    1.344002187252044678e-01,
    4.664303064346313477e-01,
    1.200614497065544128e-05,
    -2.378878416493535042e-03,
    6.179034709930419922e-01}
    })),

    B2_(scalarList({
    -7.113448381423950195e-01,
    2.472758889198303223e-01,
    3.876018822193145752e-01,
    2.375548034906387329e-01,
    -1.567474682815372944e-03,
    3.280445039272308350e-01,
    1.348470598459243774e-01
    })),

    W3_(scalarListList({
    {5.785049796104431152e-01,
    -6.650802493095397949e-01,
    -3.015841543674468994e-01,
    6.611549109220504761e-02,
    6.264905333518981934e-01},
    {-1.010075807571411133e-01,
    4.125173091888427734e-01,
    5.404242873191833496e-01,
    -2.215646207332611084e-02,
    -3.892908692359924316e-01},
    {-1.008486971259117126e-01,
    -6.076703667640686035e-01,
    -3.974271714687347412e-01,
    -5.997449532151222229e-02,
    8.935483545064926147e-02},
    {-2.227297574281692505e-01,
    1.532244682312011719e-01,
    2.081034332513809204e-01,
    -9.441027790307998657e-02,
    -7.180412411689758301e-01},
    {3.099394962191581726e-05,
    4.564229166135191917e-04,
    4.146466962993144989e-04,
    -3.556399606168270111e-04,
    2.226877841167151928e-04},
    {5.635383129119873047e-01,
    3.086688811890780926e-04,
    7.379933958873152733e-04,
    4.955981858074665070e-04,
    7.646893262863159180e-01},
    {-1.232083439826965332e+00,
    7.223671674728393555e-01,
    5.669858455657958984e-01,
    3.303100541234016418e-02,
    -6.702869664877653122e-03}
    })),

    B3_(scalarList({
    4.078924059867858887e-01,
    4.463056027889251709e-01,
    7.500135153532028198e-02,
    7.885161787271499634e-03,
    4.308234751224517822e-01
    })),

    W4_(scalarListList({{
    4.089223384857177734e+00,
    2.709713697433471680e+00,
    3.118484020233154297e+00,
    3.863466978073120117e+00,
    -3.055608272552490234e+00
    }})),

    B4_(dimensionedScalar(dimVelocity, 1.684066504240036011e-01)),

    // Input layer
    L0_(PtrList<volScalarField>(7)),
    L0_names_(List<string>(7)),
    Y0out_(PtrList<volScalarField>(7)),
    Y0out_names_(List<string>(7)),

    // Layer 1
    L1_(PtrList<volScalarField>(10)),
    L1_names_(List<string>(10)),
    Y1out_(PtrList<volScalarField>(10)),
    Y1out_names_(List<string>(10)),

    // Layer 2
    L2_(PtrList<volScalarField>(7)),
    L2_names_(List<string>(7)),
    Y2out_(PtrList<volScalarField>(7)),
    Y2out_names_(List<string>(7)),

    // Layer 3
    L3_(PtrList<volScalarField>(5)),
    L3_names_(List<string>(5)),
    Y3out_(PtrList<volScalarField>(5)),
    Y3out_names_(List<string>(5)),

    // Layer 4
    L4_
    (
        volScalarField
        (
            IOobject
            (
                "L4",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            scalar(0)
        )
    )
{
    appendInfo("\tLBV estimation method: ANN correlation");

    par_.set(0, new volScalarField("p_dimless", p_/3970000));

    par_.set
    (
        1, 
        new volScalarField
        (
            IOobject
            (
                "ER",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            ER_/7.16
        )
    );

    par_.set
    (
        2,
        new volScalarField
        (
            IOobject
            (
                "TU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("TU", dimTemperature, Zero)
        )
    );

    par_[0].dimensions().reset(dimless);

    // Input layer setup
    L0_names_[0] = "L0_0";
    L0_names_[1] = "L0_1";
    L0_names_[2] = "L0_2";
    L0_names_[3] = "L0_3";
    L0_names_[4] = "L0_4";
    L0_names_[5] = "L0_5";
    L0_names_[6] = "L0_6";

    Y0out_names_[0] = "Y0out_0";
    Y0out_names_[1] = "Y0out_1";
    Y0out_names_[2] = "Y0out_2";
    Y0out_names_[3] = "Y0out_3";
    Y0out_names_[4] = "Y0out_4";
    Y0out_names_[5] = "Y0out_5";
    Y0out_names_[6] = "Y0out_6";

    forAll(L0_names_, s)
    {
        L0_.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    L0_names_[s],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    for (int i = 0; i < 7; i++) 
    { 
        Y0out_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    Y0out_names_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    // Layer 1 setup
    L1_names_[0] = "L1_0";
    L1_names_[1] = "L1_1";
    L1_names_[2] = "L1_2";
    L1_names_[3] = "L1_3";
    L1_names_[4] = "L1_4";
    L1_names_[5] = "L1_5";
    L1_names_[6] = "L1_6";
    L1_names_[7] = "L1_7";
    L1_names_[8] = "L1_8";
    L1_names_[9] = "L1_9";

    Y1out_names_[0] = "Y1out_0";
    Y1out_names_[1] = "Y1out_1";
    Y1out_names_[2] = "Y1out_2";
    Y1out_names_[3] = "Y1out_3";
    Y1out_names_[4] = "Y1out_4";
    Y1out_names_[5] = "Y1out_5";
    Y1out_names_[6] = "Y1out_6";
    Y1out_names_[7] = "Y1out_7";
    Y1out_names_[8] = "Y1out_8";
    Y1out_names_[9] = "Y1out_9";

    forAll(L1_names_, s)
    {
        L1_.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    L1_names_[s],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    for (int i = 0; i < 10; i++) 
    { 
        Y1out_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    Y1out_names_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    // Layer 2 setup
    L2_names_[0] = "L2_0";
    L2_names_[1] = "L2_1";
    L2_names_[2] = "L2_2";
    L2_names_[3] = "L2_3";
    L2_names_[4] = "L2_4";
    L2_names_[5] = "L2_5";
    L2_names_[6] = "L2_6";

    Y2out_names_[0] = "Y2out_0";
    Y2out_names_[1] = "Y2out_1";
    Y2out_names_[2] = "Y2out_2";
    Y2out_names_[3] = "Y2out_3";
    Y2out_names_[4] = "Y2out_4";
    Y2out_names_[5] = "Y2out_5";
    Y2out_names_[6] = "Y2out_6";

    forAll(L2_names_, s)
    {
        L2_.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    L2_names_[s],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    for (int i = 0; i < 7; i++) 
    { 
        Y2out_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    Y2out_names_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    // Layer 3 setup 
    L3_names_[0] = "L3_0";
    L3_names_[1] = "L3_1";
    L3_names_[2] = "L3_2";
    L3_names_[3] = "L3_3";
    L3_names_[4] = "L3_4";

    Y3out_names_[0] = "Y3out_0";
    Y3out_names_[1] = "Y3out_1";
    Y3out_names_[2] = "Y3out_2";
    Y3out_names_[3] = "Y3out_3";
    Y3out_names_[4] = "Y3out_4";

    forAll(L3_names_, s)
    {
        L3_.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    L3_names_[s],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }

    for (int i = 0; i < 5; i++) 
    { 
        Y3out_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    Y3out_names_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                scalar(0)
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarBurningVelocityModels::ANN::~ANN()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::laminarBurningVelocityModels::ANN::correct
()
{
    if (debug_)
    {
        Info << "\t\t\tANN correct:" << endl;
        Info << "\t\t\t\tInitial average S_L: "  << average(sLaminar_).value() << endl;
    }

    // TU parameter setup
    par_[2].dimensions().reset(dimTemperature);
    par_[2] = reactionRate_.TU() / 864.0;
    par_[2].dimensions().reset(dimless);

    // Input layer calculations
    for (int i = 0; i < 7; i++) 
    {   
        L0_[i] = 0;
        for (int k = 0; k < 3; k++) 
        {
            L0_[i] += W0_[k][i] * par_[k];
        }
        Y0out_[i] = tanh(L0_[i] + B0_[i]);
    }

    // Layer 1 calculations
    for (int i = 0; i < 10; i++) 
    {   
        L1_[i] = 0;
        for (int k = 0; k < 7; k++) 
        {
            L1_[i] += W1_[k][i] * Y0out_[k];
        }
        Y1out_[i] = L1_[i] + B1_[i];
        Y1out_[i].max(0);
    }

    // Layer 2 calculations  
    for (int i = 0; i < 7; i++) 
    {   
        L2_[i] = 0;
        for (int k = 0; k < 10; k++) 
        {
            L2_[i] += W2_[k][i] * Y1out_[k];
        }
        Y2out_[i] = L2_[i] + B2_[i];
        Y2out_[i].max(0);
    }

    // Layer 3 calculations
    for (int i = 0; i < 5; i++) 
    {   
        L3_[i] = 0;
        for (int k = 0; k < 7; k++) 
        {
            L3_[i] += W3_[k][i] * Y2out_[k];
        }
        Y3out_[i] = L3_[i] + B3_[i];
        Y3out_[i].max(0);
    }
        
    // Layer 4 setup
    L4_.dimensions().reset(dimless);
    L4_ = 0;
        
    // Layer 4 calculations 
    for (int k = 0; k < 5; k++)
    {   
        L4_ += W4_[0][k] * Y3out_[k];
    }

    L4_.dimensions().reset(dimVelocity);
    sLaminar_ = L4_ + B4_;
    sLaminar_.max(0);

    if (debug_)
    {
        Info << "\t\t\t\tObtained average S_L: "  << average(sLaminar_).value() << endl;
        Info << "\t\t\t\tANN correct finished" << endl;
    }

     
}
// ************************************************************************* //
