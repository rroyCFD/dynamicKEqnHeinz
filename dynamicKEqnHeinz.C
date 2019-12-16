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
volScalarField dynamicKEqnHeinz<BasicTurbulenceModel>::Ck
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volSymmTensorField LL
    (
        simpleFilter_(dev(filter_(sqr(this->U_)) - (sqr(filter_(this->U_)))))
    );

    const volSymmTensorField MM
    (
        simpleFilter_(-2.0*this->delta()*sqrt(KK)*filter_(D))
    );

    const volScalarField Ck
    (
        simpleFilter_(0.5*(LL && MM))
       /(
            simpleFilter_(magSqr(MM))
          + dimensionedScalar("small", sqr(MM.dimensions()), vSmall)
        )
    );

    tmp<volScalarField> tfld = 0.5*(mag(Ck) + Ck);
    return tfld();
}


template<class BasicTurbulenceModel>
volScalarField dynamicKEqnHeinz<BasicTurbulenceModel>::Ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volScalarField Ce
    (
        simpleFilter_(this->nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D))))
       /simpleFilter_(pow(KK, 1.5)/(2.0*this->delta()))
    );

    tmp<volScalarField> tfld = 0.5*(mag(Ce) + Ce);
    return tfld();
}


template<class BasicTurbulenceModel>
volScalarField dynamicKEqnHeinz<BasicTurbulenceModel>::Ce() const
{
    const volSymmTensorField D(dev(symm(fvc::grad(this->U_))));

    volScalarField KK
    (
        0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_)))
    );
    KK.max(dimensionedScalar("small", KK.dimensions(), small));

    return Ce(D, KK);
}


template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctNut
(
    const volSymmTensorField& D,
    const volScalarField& KK
)
{
    this->nut_ = Ck(D, KK)*sqrt(k_)*this->delta();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctNut()
{
    const volScalarField KK
    (
        0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_)))
    );

    correctNut(symm(fvc::grad(this->U_)), KK);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> dynamicKEqnHeinz<BasicTurbulenceModel>::kSource() const
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
    filterRatio_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "filterRatio",
            this->coeffDict_,
            2.0
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
    dynamic_
    (
        Switch::lookupOrAddToDict
        (
            "dynamic",
            this->coeffDict_,
            true
        )
     ),
    damping_
    (
        Switch::lookupOrAddToDict
        (
            "damping",
            this->coeffDict_,
            true
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
    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, this->kMin_);

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
            Ce()*k()*sqrt(k())/this->delta()
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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volSymmTensorField D(dev(symm(tgradU())));
    const volScalarField G(this->GName(), 2.0*nut*(tgradU() && D));
    tgradU.clear();

    volScalarField KK(0.5*(filter_(magSqr(U)) - magSqr(filter_(U))));
    KK.max(dimensionedScalar("small", KK.dimensions(), small));

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
    ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(Ce(D, KK)*alpha*rho*sqrt(k_)/this->delta(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(D, KK);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
