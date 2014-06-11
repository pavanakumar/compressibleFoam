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

/// Loop over internal faces and get face flux from L and R cell center
label leftCell , rightCell;
forAll( mesh.owner() , iface ) {
  /// Get the left and right cell index
  leftCell = mesh.owner()[iface];
  rightCell = mesh.neighbour()[iface]; 
  /// Approximate Riemann solver at interface 
  scalar lambda = (*fluxSolver)( &rho[leftCell]  , &U[leftCell]  , &p[leftCell] , /// Left
       &rho[rightCell] , &U[rightCell] , &p[rightCell] ,          /// Right
       &massFlux[iface] , &momFlux[iface] , &energyFlux[iface] ,
       &nf[iface] );
  /// Multiply with face area to get face flux
  massFlux[iface] *= mesh.magSf()[iface];
  momFlux[iface] *= mesh.magSf()[iface];
  energyFlux[iface] *= mesh.magSf()[iface];
  localDt[ leftCell ]  += lambda * mesh.magSf()[iface];
  localDt[ rightCell ] += lambda * mesh.magSf()[iface];
}

