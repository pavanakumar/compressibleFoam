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

/// Default flux solver is roe scheme 
word fluxScheme("roe");
fluxSolver = &roe;

if (mesh.schemesDict().found("fluxScheme"))
{
  fluxScheme = word(mesh.schemesDict().lookup("fluxScheme"));
  if ((fluxScheme == "roe") ){
    Info<< "fluxScheme: " << fluxScheme << endl;
    fluxSolver = &roe;
  }
  else if( (fluxScheme == "van_leer") ){
    Info<< "fluxScheme: " << fluxScheme << endl;
    fluxSolver = &van_leer;
  }
  else{
    FatalErrorIn("EulerSolver::readFluxScheme") 
      << "fluxScheme: " << fluxScheme
      << " is not a valid choice. "
      << "Options are: roe,van_leer"
      << abort(FatalError);
  }
}


