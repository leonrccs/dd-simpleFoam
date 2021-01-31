/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "myIterativeTurbulenceModel.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void myIterativeTurbulenceModel<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> myIterativeTurbulenceModel<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> myIterativeTurbulenceModel<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
myIterativeTurbulenceModel<BasicTurbulenceModel>::myIterativeTurbulenceModel
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    gamma_start_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma_start",
            this->coeffDict_,
            0.0
        )
    ),
    gamma_max_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma_max",
            this->coeffDict_,
            0.1
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    b_ml
    (
        IOobject
        (
            IOobject::groupName("b_ml", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    a_ml
    (
        IOobject
        (
            "a_ml",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("a_ml", dimensionSet(0, 2, -2, 0, 0, 0, 0), symmTensor(.0, .0, .0, .0, .0, .0)) //tensor(0, 0, 0, 0, 0, 0, 0 ,0 ,0))
    ),
    a_rans
    (
        IOobject
        (
            "a_rans",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("a_rans", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor(.0, .0, .0, .0, .0, .0, .0, .0, .0)) //tensor(0, 0, 0, 0, 0, 0, 0 ,0 ,0))
    ),
    a_mix
    (
        IOobject
        (
            "a_mix",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("a_mix", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor(.0, .0, .0, .0, .0, .0, .0, .0, .0)) //tensor(0, 0, 0, 0, 0, 0, 0 ,0 ,0))
    ),
    gamma_mix_it_
    (
        dimensionedScalar("gamma_mix_it", dimless, 0.0)
    ),
    ident_
    (
        IOobject
        (
            "ident_",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
    ),
    gamma_mix
    (
        IOobject
        (
            "gamma_mix",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("gamma_mix", dimless, 0.0)
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool myIterativeTurbulenceModel<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        gamma_max_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}

template<class BasicTurbulenceModel>
void myIterativeTurbulenceModel<BasicTurbulenceModel>::updateGamma()
{
    const float thisTime = this->runTime_.value();
    const float lastTime = readScalar(this->runTime_.controlDict().lookup("endTime"));
    gamma_mix_it_ = gamma_start_ + (gamma_max_ - gamma_start_) * (thisTime/lastTime);
    Info << "Current time:       " << thisTime << endl;
    Info << "End time:           " << readScalar(this->runTime_.controlDict().lookup("endTime")) << endl;
    Info << "Mixing coefficient: " << thisTime/lastTime << endl;
    Info << "Updating mixing coefficient. Current mixing level:  " << gamma_mix_it_ << endl;
}

template<class BasicTurbulenceModel>
void myIterativeTurbulenceModel<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // upate mixing coefficient
    updateGamma();

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        this->type() + ":GbyNu",
        tgradU().v() && dev(twoSymm(tgradU().v()))
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);

    const volScalarField::Internal G_ml(
        this->GName(),
        //tgradU().v() && (this->gamma_mix.v()*this->a_ml.v())
        tgradU().v() && this->a_ml.v()
    );

    const volScalarField::Internal G_mixed(
        this->GName(),
        (1 - gamma_mix_it_)*G - gamma_mix_it_*G_ml
    );

    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*GbyNu*Cmu_*k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G_mixed //standard
        //(1 - this->gamma_mix)*alpha()*rho()*G
      //- this->gamma_mix*alpha()*rho()*G_ml
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}



template<class BasicTurbulenceModel>
tmp<fvVectorMatrix>myIterativeTurbulenceModel<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    // upate mixing coefficient
    const float thisTime = this->runTime_.value();
    const float startTime = readScalar(this->runTime_.controlDict().lookup("startTime"));
    const float endTime = readScalar(this->runTime_.controlDict().lookup("endTime"));
    gamma_mix_it_ = gamma_start_ + (gamma_max_ - gamma_start_) * ((thisTime - startTime)/(endTime - startTime));
    Info << "Updating mixing coefficient. Current mixing level:  " << gamma_mix_it_.value() << endl;

    // calculate anisotropy tensor
    this->a_ml = 2*this->alpha_*this->rho_*this->k()*(this->b_ml);
    this->a_rans = -1*this->alpha_*this->rho_*this->nut_*(ident_ & twoSymm(fvc::grad(U))); // & 
    this->a_mix = (1 - gamma_mix_it_)*this->a_rans + gamma_mix_it_*this->a_ml;

    // Mixing field
    // Defined from https://doi.org/10.2514/1.2094
    float lh = 1.0; // Characteristic length
    volScalarField d = wallDist(this->mesh_).y()/dimensionedScalar("lh", dimensionSet(0, 1, 0, 0, 0, 0, 0), lh);
    volScalarField lambdb0 = sqrt(this->k()/this->epsilon());
    // float lambdb0 = 1.0;
    lambdb0.dimensions().reset(dimless);

    // TODO: read scaling parameters from tubulence model parameters
    dimensionedScalar alpha1 = dimensionedScalar("alpha1", dimless, 1.0);
    dimensionedScalar alpha2 = dimensionedScalar("alpha2", dimless, 0.05);
    dimensionedScalar small0 = dimensionedScalar("small", dimless, 1e-10);

    // mixing coefficient
    // Geneva's implementation:
    this->gamma_mix = min(alpha1, pow(tanh((d/(alpha2*lambdb0 + small0))), 2));
    //Info << "Mixing coefficient at cell 1000: " << this->gamma_mix.internalField()[1000] << endl;
    //Info << "Max mixing coefficient: " << max(this->gamma_mix.internalField()) << endl;
    // Kaandorp's implementation:

    return
    (
    // Original k-epsilon formulation:
    //  - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
    //  - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    //  Split up into viscous and eddy viscosity formulation:
    //  - fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)
    //  - fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
    //  - fvm::laplacian(this->alpha_*this->rho_*this->nut_, U)
    //  - fvc::div((this->alpha_*this->rho_*this->nut_)*dev2(T(fvc::grad(U))))
    // Including data driven anisotropy tensor and genenevas gamma mix
    //   - fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)
    //   - fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
    //   - fvm::laplacian(this->alpha_*this->rho_*this->nut_*(1-this->gamma_mix), U)
    //   - fvc::div((this->alpha_*this->rho_*this->nut_)
    //   *dev2(T(fvc::grad(U)))*(1-this->gamma_mix))
    //   + fvc::div(this->gamma_mix*this->a_ml)
    // Inlcuding data driven anisotropy tensor and kaandorps approach
      - fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)
      - fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nut_*(1 - this->gamma_mix * gamma_mix_it_), U)
      - fvc::div((this->alpha_*this->rho_*this->nut_)
      *dev2(T(fvc::grad(U)))*(1 - this->gamma_mix * gamma_mix_it_))
      + fvc::div(this->gamma_mix * gamma_mix_it_*this->a_ml)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
