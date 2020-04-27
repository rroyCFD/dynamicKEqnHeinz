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

Class
     Foam::LESModels::dynamicKEqnHeinz

Description
    Dynamic one equation eddy-viscosity model

    Eddy viscosity SGS model using a modeled balance equation to simulate
    the behaviour of k in which a dynamic procedure is applied to evaluate the
    coefficients.

    Reference:
    \verbatim
        Kim, W and Menon, S. (1995).
        A new dynamic one-equation subgrid-scale model for
        large eddy simulation.
        In 33rd Aerospace Sciences Meeting and Exhibit, Reno, NV, 1995.
    \endverbatim

    There are no default model coefficients but the filter used for KK must be
    supplied, e.g.
    \verbatim
        dynamicKEqnHeinzCoeffs
        {
            filter simple;
        }
    \endverbatim

SourceFiles
    dynamicKEqnHeinz.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicKEqnHeinz_H
#define dynamicKEqnHeinz_H

#include "LESeddyViscosity.H"
#include "simpleFilter.H"
#include "bound.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynamicKEqnHeinz Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class dynamicKEqnHeinz
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        dynamicKEqnHeinz(const dynamicKEqnHeinz&);
        void operator=(const dynamicKEqnHeinz&);

protected:

    // Protected data

        // Model constants
            dimensionedScalar Ckmin_;
            dimensionedScalar Ckmax_;
            dimensionedScalar Ce_;

            dimensionedScalar filterRatio_;

            Switch calcMoreFields_;

        // Fields

            volScalarField k_;
            volScalarField Ckd_;
            volSymmTensorField B_;

        // Filters

            simpleFilter simpleFilter_;
            autoPtr<LESfilter> filterPtr_;
            LESfilter& filter_;


    // Protected Member Functions

        //- Update dynamic Ck coefficient field
        virtual void updateCkd(tmp<volTensorField>&);

        //- Update sub-grid eddy-viscosity
        virtual void correctNut();

        //- Update SGS stress (write on write only)
        virtual void correctB();
        virtual void correctB(tmp<volTensorField>& );

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("dynamicKEqnHeinz");


    // Constructors

        //- Construct from components
        dynamicKEqnHeinz
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~dynamicKEqnHeinz()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return sub-grid scale stress
        virtual tmp<volSymmTensorField> B() const
        {
            return B_;
        }

        // Used in the RHS of momentum equation.
        tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const
        {
            return
            (
                - fvm::laplacian(this->nuEff_, U)
                - fvc::div(this->nuEff_*dev(T(fvc::grad(U))))
            );
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", this->nut_ + this->nu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicKEqnHeinz.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
