/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "dynamicKEqnHeinz.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctNut()
{
    this->nut_= Ckd_* this->delta() * sqrt(k_);

    // If non-dynamic, apply damping function if required
    if(!dynamic_ && damping_)
    {
        // damping function
        volScalarField Ret = this->nut_ / this->nu();
        volScalarField fmu_ =
          0.09 + (0.91 + 1./pow(Ret+SMALL,3)) * (1.0 - exp(-pow(Ret/25.,2.75)));
        this ->nut_ *= fmu_;
    }

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctB()
{
    B_= ((2.0/3.0)*I)*k_ - 2.0*this->nut_*symm(dev(fvc::grad(this->U_))) + N_;
    B_.correctBoundaryConditions();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicKEqnHeinz<BasicTurbulenceModel>::dynamicKEqnHeinz
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
    LESeddyViscosity<BasicTurbulenceModel>
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
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.1
        )
    ),
    Ckmin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ckmin",
            this->coeffDict_,
            -0.05
        )
    ),
    Ckmax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ckmax",
            this->coeffDict_,
            0.5
        )
    ),
    Cnmin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnmin",
            this->coeffDict_,
            0.
        )
    ),
    Cnmax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnmax",
            this->coeffDict_,
            5.0
        )
    ),
    Ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ce",
            this->coeffDict_,
            1.048
        )
    ),
    filterRatio_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "filterRatio",
            this->coeffDict_,
            2.0
        )
    ),
    dynamic_
    (
        Switch::lookupOrAddToDict
        (
            "dynamic",
            this->coeffDict_,
            true
        )
    ),
    nonLinear_
    (
        Switch::lookupOrAddToDict
        (
            "nonLinear",
            this->coeffDict_,
            false
        )
    ),
    damping_
    (
        Switch::lookupOrAddToDict
        (
            "damping",
            this->coeffDict_,
            false
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Ckd_
    (
        IOobject
        (
            IOobject::groupName("Ckd", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Ckd", dimless, SMALL)
    ),
    Cnd_
    (
        IOobject
        (
            IOobject::groupName("Cnd", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Cnd", dimless, SMALL)
    ),
    N_
    (
        IOobject
        (
            IOobject::groupName("N", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    B_
    (
        IOobject
        (
            IOobject::groupName("B", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor
            ("B", dimensionSet(0, 2, -2, 0, 0, 0, 0), symmTensor::zero)
    ),
    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, this->kMin_);

    // update SGS stress B
    correctB();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool dynamicKEqnHeinz<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {

        Ckmin_.readIfPresent(this->coeffDict());
        Ckmax_.readIfPresent(this->coeffDict());
        Cnmin_.readIfPresent(this->coeffDict());
        Cnmax_.readIfPresent(this->coeffDict());
        Ce_.readIfPresent(this->coeffDict());
        filterRatio_.readIfPresent(this->coeffDict());

        nonLinear_.readIfPresent("nonLinear", this->coeffDict());
        dynamic_.readIfPresent("dynamic", this->coeffDict());
        damping_.readIfPresent("damping", this->coeffDict());

        filter_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> dynamicKEqnHeinz<BasicTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Ce_*k()*sqrt(k())/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correct()
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
    // volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));

    /*
    S_ij-d and Omega_ij:

    a) Important!!!! OpenFOAM defines grad(U) = dU_j/dxi -> hence we have to
    use grad(U)^T (transpose) to make it match our Cartesian tensor notation!
    b) gradU.T() should be trace-less due to continuity equation but deviates
    slightly due to numerics
    */

    const volSymmTensorField Sijd (dev(symm(tgradU())));
    const volTensorField     Rotij(skew(tgradU().T()));

    // filtered S_ij-d and magnitude
    // (careful -> use Stefans deifinition abs(L) = sqrt(2 Lij Lji)
    const volSymmTensorField SijdF  = dev(filter_(Sijd)); // this is D
    const volScalarField magSf = sqrt(2.) * mag(SijdF);

    // Leonard stress --------------------------------------------------------//
    volSymmTensorField Lijd = (filter_(sqr(U)) - (sqr(filter_(U))));
    // Test-filter kinetic energy
    const volScalarField ktest = 0.5 * tr(Lijd); // this is KK
    // Deviatoric part of Leonard stress and it's magnitude
    Lijd = dev(Lijd);
    const volScalarField magLd = sqrt(2.) * mag(Lijd);

    const volScalarField G(this->GName(), 0.5*tr(-twoSymm(B() & tgradU())));
    tgradU.clear();


    if (dynamic_)
    {
        // Correlation coeffcients
        dimensionedScalar small1
        (
            "small1",
            dimensionSet(0, 2, -3, 0, 0, 0, 0),
            SMALL
        );

        // Test-filter width
        volScalarField deltaT = filterRatio_*this->delta();

        // Non-Linear Model
        if (nonLinear_)
        {
            // Tensors for dynamic procedure
            const volSymmTensorField nij =
            (
                sqr(deltaT) *
                (
                    twoSymm(SijdF & filter_(Rotij))
                    - twoSymm(SijdF & SijdF)
                    + 2./3.*I*(SijdF && SijdF)
                )
            );
            const volScalarField magn = sqrt(2.) * mag(nij);


            // Correlation coeffcients
            dimensionedScalar small2
            (
                "small2",
                dimensionSet(0, 4, -4, 0, 0, 0, 0),
                SMALL
            );

            const volScalarField rSN = (SijdF && nij)  / (0.5*magSf*magn + small1);
            const volScalarField rLN = (Lijd  && nij)  / (0.5*magLd*magn + small2);
            const volScalarField rLS = (Lijd  && SijdF)/ (0.5*magLd*magSf + small1);

            // Calculate the dynamic coeffcients
            Ckd_.primitiveFieldRef() =
            (
                (rSN*rLN - rLS) * magLd /
                    (
                        (2.*(deltaT*sqrt(max(ktest, this->kMin_*0.)))
                        *magSf*(1. - sqr(rSN)))
                    + this->kMin_
                    )
            );

            Cnd_.primitiveFieldRef() =
                (rSN*rLS - rLN) * magLd / ((1. - sqr(rSN))*magn + this->kMin_);

            // Clipping of Ck
            Ckd_.max(Ckmin_);
            Ckd_.min(Ckmax_);

            // Clipping of Cn
            Cnd_.max(Cnmin_);
            Cnd_.min(Cnmax_);
        }
        else
        {
            const volScalarField rLS =
                (Lijd  && SijdF)/ (0.5*magLd*magSf + small1);

            // Calculate the dynamic coeffcients
            Ckd_.primitiveFieldRef() =
            (
                (-rLS) * magLd /
                ((2.*(deltaT*sqrt(max(ktest,this->kMin_*0.)))*magSf)
                + this->kMin_)
            );

            Cnd_.primitiveFieldRef() = 0.* Ckd_;

            // Clipping of Ck
            Ckd_.max(Ckmin_);
            Ckd_.min(Ckmax_);
        }
    }
    else
    {
        if (nonLinear_)
        {
            Ckd_ = Ck_;
            Cnd_ = 3.*sqr(Ck_);
        }
        else
        {
            Ckd_ = Ck_;
            Cnd_ = 0.;
        }
    }

    Info<< "Constant: Ck:"<< max(Ckd_).value()<< tab<< min(Ckd_).value()<< endl;
    Info<< "Constant: Cn:"<< max(Cnd_).value()<< tab<< min(Cnd_).value()<< endl;


    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
    ==
        alpha*rho*G
      - fvm::Sp(Ce_*alpha*rho*sqrt(k_)/this->delta(), k_) // fix Ce()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    // update SGS viscosity
    correctNut();

    // update NonLinear Stress
    N_ = - Cnd_ * sqr(this->delta()) *
        (twoSymm(Sijd & Rotij) - twoSymm(Sijd & Sijd) + 2./3.*I*(Sijd && Sijd));
    N_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
