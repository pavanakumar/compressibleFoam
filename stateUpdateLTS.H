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
  rhoResidMax = 0.0;
  scalar lambda_ds_max = 0.0,
         lambda_ds_min = std::numeric_limits<scalar>::max();
  forAll( rho , icell ){
    if( fabs( localDt[icell] ) > lambda_ds_max )
      lambda_ds_max = fabs( localDt[icell] );
    if( fabs( localDt[icell] ) < lambda_ds_min )
      lambda_ds_min = fabs( localDt[icell] );
  } 
  /// Convert to conservative
  forAll( rho , icell ){
    convertConservative( &rho[icell] , &U[icell] , &p[icell] );
    localDt[icell] = CFL / localDt[icell] ;
    rho[icell] = rho[icell] - localDt[icell] * massResidue[icell];
    U[icell]   = U[icell]   - localDt[icell] * momResidue[icell];
    p[icell]   = p[icell]   - localDt[icell] * energyResidue[icell];
    /// Calcuate maximum residue
    if( fabs( localDt[icell] * massResidue[icell] ) > rhoResidMax ) 
      rhoResidMax = fabs( localDt[icell] * massResidue[icell] );
    /// Convert to primitive
    convertPrimitive( &rho[icell] , &U[icell] , &p[icell] );
    localDt[icell] = 0.0;
  }

  Info << "Lambda min/max = " << lambda_ds_max << "  " << lambda_ds_min << "  ";

  if ( Pstream::parRun() )
    reduce( rhoResidMax , maxOp<scalar>() );
   
