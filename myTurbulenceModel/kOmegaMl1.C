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

#include "kOmegaMl1.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaMl1<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaMl1<BasicTurbulenceModel>::kOmegaMl1
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
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.072
        )
    ),
    gamma_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
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
    end_ramp_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "end_ramp",
            this->coeffDict_,
            0.5
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
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
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
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaMl1<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kOmegaMl1<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

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

    // get current mixing coefficient
    dimensioned<scalar> thisTime = this->runTime_.value();
    dimensioned<scalar> startTime = readScalar(this->runTime_.controlDict().lookup("startTime"));
    dimensioned<scalar> endTime = readScalar(this->runTime_.controlDict().lookup("endTime"));
    dimensioned<scalar> rampTime = end_ramp_*(endTime - startTime) + startTime;
    dimensioned<scalar> gamma_mix_it_ = gamma_start_ + (gamma_max_ - gamma_start_) * ((thisTime - startTime)/(rampTime - startTime));

    if (thisTime > rampTime)
    {
        gamma_mix_it_ = gamma_max_;
    }

    Info << "Current mixing coefficient:" << gamma_mix_it_.value() << endl;

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        tgradU().v() && dev(twoSymm(tgradU().v()))
    );
    // const volScalarField::Internal G(this->GName(), nut()*GbyNu);

    const volScalarField::Internal G
    (
        this->GName(),
        // nut *( tgradU () && dev ( twoSymm ( tgradU () )))
        gamma_mix_it_*(-2)*k_*((b_ml+(1.0/3.0)*I)&&tgradU())
        + (1 - gamma_mix_it_)*nut*(tgradU() && dev(twoSymm(tgradU())))
    );

    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*G*omega_/k_  //gamma_*alpha()*rho()*GbyNu
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(beta_*alpha()*rho()*omega_(), omega_)
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(Cmu_*alpha()*rho()*omega_(), k_)
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
tmp<fvVectorMatrix>kOmegaMl1<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    // upate mixing coefficient
    dimensioned<scalar> thisTime = this->runTime_.value();
    dimensioned<scalar> startTime = readScalar(this->runTime_.controlDict().lookup("startTime"));
    dimensioned<scalar> endTime = readScalar(this->runTime_.controlDict().lookup("endTime"));
    dimensioned<scalar> rampTime = end_ramp_*(endTime - startTime) + startTime;
    dimensioned<scalar> gamma_mix_it_ = gamma_start_ + (gamma_max_ - gamma_start_) * ((thisTime - startTime)/(rampTime - startTime));

    if (thisTime > rampTime)
    {
        gamma_mix_it_ = gamma_max_;
    }

    Info << "Computing Reynolds stresses from NN prediction" << endl;

    return
    (
    // Original eddy viscosity formulation:
    //  - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
    //  - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    // Data driven b incorporated
     - (1-gamma_mix_it_)*fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
     - (1-gamma_mix_it_)*fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
     + gamma_mix_it_*fvc::div(2*this->k_ * ((b_ml) + (1.0/3.0)*I))
     - gamma_mix_it_*fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
     - gamma_mix_it_*fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
