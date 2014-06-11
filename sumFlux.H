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

   /// Clear all residues
   massResidue = dimensionedScalar( "", massResidue.dimensions() , 0.0 );
   momResidue = dimensionedVector( "" , momResidue.dimensions() , vector( 0.0 , 0.0 , 0.0 ) );
   energyResidue = dimensionedScalar( "" , energyResidue.dimensions() ,  0.0 );
   /// Loop over each face and add the flux to left and subtract from the right
   forAll( mesh.owner() , iface ) {
     /// Store left and right cell reference
     leftCell = mesh.owner()[iface];
     rightCell = mesh.neighbour()[iface];
     /// Note that the normal vector to a face will point
     /// from the owner cell to the neighbour cell (L to R)
     /// Add to left cell
     massResidue[leftCell]    += massFlux[iface];
     momResidue[leftCell]     += momFlux[iface];
     energyResidue[leftCell]  += energyFlux[iface];
     /// Subtract from right cell
     massResidue[rightCell]   -= massFlux[iface];
     momResidue[rightCell]    -= momFlux[iface];
     energyResidue[rightCell] -= energyFlux[iface];
   }
