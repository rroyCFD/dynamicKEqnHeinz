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
void dynamicKEqnHeinz<BasicTurbulenceModel>::updateCkd(
                                                    tmp<volTensorField>& tgradU)
{
    tmp<volVectorField> tUf = filter_(this->U_);
    volVectorField& Uf = tUf.ref();

    // filtered S_ij-d and magnitude
    // (careful -> use Stefans deifinition abs(L) = sqrt(2 Lij Lji)
    const volSymmTensorField SijdF(dev(symm(fvc::grad(Uf))));
    const volScalarField magSdf = sqrt(2.) * mag(SijdF);

    // Leonard stress --------------------------------------------------------//
    volSymmTensorField Lijd = (filter_(sqr(this->U_)) - (sqr(Uf)));
    if(calcMoreFields_ && this->runTime_.outputTime())
    {
        Uf.rename("U-SGSFilter");
        Uf.write();
    }
    tUf.clear();

    // Test-filter kinetic energy
    const volScalarField kTest("kTest", 0.5 * tr(Lijd));
    if(calcMoreFields_ && this->runTime_.outputTime())
    {
        kTest .write();
    }

    // Deviatoric part of Leonard stress and it's magnitude
    Lijd = dev(Lijd);
    const volScalarField magLd = sqrt(2.) * mag(Lijd);


    // Correlation coeffcients
    dimensionedScalar small1
    (
        "small1",
        dimensionSet(0, 2, -3, 0, 0, 0, 0),
        SMALL
    );

    // Test-filter width
    volScalarField deltaT = filterRatio_*this->delta();

    const volScalarField rLS =
        (Lijd  && SijdF)/ (0.5*magLd*magSdf + small1);

    // Calculate the dynamic coeffcients
    Ckd_.primitiveFieldRef() =
    (
        (-rLS) * magLd /
        ((2.*(deltaT*sqrt(max(kTest ,this->kMin_*0.)))*magSdf)
        + this->kMin_)
    );

    // Clipping of Ck
    Ckd_.max(Ckmin_);
    Ckd_.min(Ckmax_);

    Info << "Ck min:" << gMin(Ckd_) << "\tmax: " << gMax(Ckd_) << endl;
}


template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctNut()
{
    this->nut_= Ckd_* this->delta() * sqrt(k_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctB()
{
    Info << "Initializing SGS stress B" << endl;
    B_= ((2.0/3.0)*I)*k_ - 2.0*this->nut_*symm(dev(fvc::grad(this->U_)));
    B_.correctBoundaryConditions();
}

template<class BasicTurbulenceModel>
void dynamicKEqnHeinz<BasicTurbulenceModel>::correctB
(tmp<volTensorField>& tgradU)
{
    Info << "Correcting SGS stress B" << endl;
    B_= ((2.0/3.0)*I)*k_ - 2.0*this->nut_*symm(dev(tgradU));
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
    calcMoreFields_
    (
        Switch::lookupOrAddToDict
        (
            "calcMoreFields",
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
        dimensionedSymmTensor ("B", k_.dimensions(), symmTensor::zero)
    ),
    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, this->kMin_);

    // update SGS stress tensor B
    if(calcMoreFields_)
    {
        correctB();
    }

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
        Ce_.readIfPresent(this->coeffDict());
        filterRatio_.readIfPresent(this->coeffDict());
        calcMoreFields_.readIfPresent("calcMoreFields", this->coeffDict());

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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));

    // update dynamic Ck coefficient field
    updateCkd(tgradU);

    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));
    // tgradU.clear();

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
    ==
        alpha*rho*G
      - fvm::Sp(Ce_*alpha*rho*sqrt(k_)/this->delta(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // update SGS viscosity
    correctNut();

    // update SGS stress tensor
    if(calcMoreFields_)
    {
        correctB(tgradU);
    }

    tgradU.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
