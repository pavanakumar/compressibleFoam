/*---------------------------------------------------------------------------*\
    eulerSolver - A simple parallel first order Gas-dynamics solver
                  based on the OpenFOAM library
    Copyright (C) 2012, Pavanakumar Mohanamuraly

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/// Derived variables
scalar rho_inf = 1.0, p_inf = 1.0 / ( gama * M_inf );
scalar rhoResidMax = 0.0;
vector u_inf ( 1.0, 0.0, 0.0 );

//////////////////////
/// Fluid Density  ///
//////////////////////
volScalarField rho
(
  IOobject
  (
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar( "rho", dimless, rho_inf )
);
//////////////////////
/// Fluid Velocity ///
//////////////////////
volVectorField U
(
  IOobject
  (
    "U",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector( "U", dimless, u_inf )
);
//////////////////////
/// Fluid Pressure ///
//////////////////////
volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar( "p", dimless, p_inf )
);
//////////////////////////////////////////////
/// Face Momentum Flux ( \rhouu \cdot ds ) ///
//////////////////////////////////////////////
surfaceVectorField momFlux
(
    IOobject
    (
        "momFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector( "momFlux", dimless, vector( 0.0 , 0.0 , 0.0 ) )
);
//////////////////////////////////////////
/// Face Mass Flux ( \rho u \cdot ds ) ///
//////////////////////////////////////////
surfaceScalarField massFlux
(
    IOobject
    (
        "massFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "massFlux", dimless, 0.0 )
);
////////////////////////
/// Face Energy Flux ///
////////////////////////
surfaceScalarField energyFlux
(
    IOobject
    (
        "energyFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "energyFlux", dimless, 0.0 )
);
//////////////////////
///  Mass residue  ///
//////////////////////
volScalarField massResidue
(
    IOobject
    (
        "massResidue",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "massResidue", dimless, 0.0)
);
////////////////////////
/// Momentum residue ///
////////////////////////
volVectorField momResidue
(
    IOobject
    (
        "momResidue",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector( "momResidue", dimless, vector( 0.0 , 0.0 , 0.0 ) )
);
//////////////////////
/// Energy residue ///
//////////////////////
volScalarField energyResidue
(
    IOobject
    (
        "energyResidue",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "energyResidue", dimless, 0.0 )
);
///////////////////////
/// Local Time Step ///
///////////////////////
volScalarField localDt
(
    IOobject
    (
        "localDt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "localDt", dimless , 0.0 )
);

// Unit face normals
surfaceVectorField nf = mesh.Sf() / mesh.magSf();
 
