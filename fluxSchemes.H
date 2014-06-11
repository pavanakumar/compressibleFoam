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

/// Constants
scalar gama = 1.4;
scalar gbygm1 = gama / ( gama - 1.0 ) , gm1 = gama - 1.0;

/// The normal flux to a face given the state vector 
/// and normal at a face
scalar normalFlux( scalar *rho , vector *u , scalar *p ,
                 vector *normal , scalar f[5] ){
  scalar ret = (*u) & (*normal);
  f[0] = (*rho) * ret;
  for( int i = 0 ; i < 3 ; ++i )
    f[i+1] = f[0] * (*u)[i] + (*p) * (*normal)[i];
  f[4] = f[0] * ( gbygm1 * (*p) / (*rho) + 0.5 * ( (*u) & (*u) ) );
  return std::fabs( ret ) + std::sqrt( gama * (*p) / (*rho) );
}

//
void convertConservative( scalar *rho , vector *u , scalar *p ){
  // \rho e_t = p / ( \gamma - 1 ) + 1/2 * \rho * u^2
  (*p) = (*p) / gm1 + 0.5 * (*rho) * ( (*u) & (*u) );
  // (\rho u,\rho v,\rho w) = (u,v,w) * \rho
  (*u) *= (*rho);
}

//
void convertPrimitive( scalar *mass , vector *mom , scalar *energy ){
  (*mom) /= (*mass);
  // p = (\gamma - 1) (\rho e_t - 1/2 \rho u^2 )
  (*energy) = gm1 * ( (*energy) - 0.5 * (*mass) * ( (*mom) & (*mom) ) );
}


/// Riemann solver that takes in density,velocity and 
/// pressure and return the flux remember that it should
/// return the scalar value of ( |u_n| + a ) at the interface
scalar (*fluxSolver)( scalar *rhoL , vector *uL , scalar *pL , 
          scalar *rhoR , vector *uR , scalar *pR ,
          scalar *mass , vector *mom , scalar *energy ,
          vector *normal ) = NULL;

/// Try implementing your own riemann flux scheme
/// This is Roe's Approximate Riemann solver
scalar roe( scalar *rhoL , vector *uL , scalar *pL ,
          scalar *rhoR , vector *uR , scalar *pR ,
          scalar *mass , vector *mom , scalar *energy ,
          vector *normal ){
  scalar h0L, h0R, Rfac, rhoRoe, aRoe,
         h0Roe, drho, dp, VnL, VnR, VnRoe, qSqRoe,
         dV, drhoTerm , uL2 , uR2;
  scalar fL[5], fR[5], df1[5], df2[5], df3[5] , dot_prod;
  vector dvel , velRoe;
  //  Left and Right Flow values
  uL2 = (*uL) & (*uL);
  uR2 = (*uR) & (*uR);
  h0L = gbygm1 * (*pL) / (*rhoL) + 0.50 * ( uL2 );
  h0R = gbygm1 * (*pR) / (*rhoR) + 0.50 * ( uR2 );
  VnL = (*uL) & (*normal);
  VnR = (*uR) & (*normal);
  //  Roe Averaged States
  Rfac = std::sqrt( (*rhoR) / (*rhoL) );
  rhoRoe = Rfac * (*rhoL);
  velRoe = ( (*uL) + Rfac * (*uR) ) / ( 1.0 + Rfac );
  h0Roe = (h0L + Rfac * h0R) / ( 1.0 + Rfac );
  qSqRoe = velRoe & velRoe ;
  aRoe = std::sqrt( gm1 * ( h0Roe - 0.50 * qSqRoe ) );
  VnRoe = (VnL + Rfac * VnR) / (1.0 + Rfac);
  //  Differential values for characteristics
  drho = (*rhoR) - (*rhoL);
  dvel = (*uR) - (*uL);
  dp   = (*pR) - (*pL);
  dV = dvel & (*normal);
  normalFlux( rhoL , uL , pL , normal , fL );
  normalFlux( rhoR , uR , pR , normal , fR );
  // Form Inteface roe averaged fluxes df1 to df3
  // Calculate df1
  drhoTerm = drho - dp / ( aRoe * aRoe );
  df1[0]   = std::abs(VnRoe) * drhoTerm;
  dot_prod = 0.0;
  for( int i = 0 ; i < 3 ; ++i ){
    df1[i+1] = std::abs(VnRoe) * ( velRoe[i] * drhoTerm +
               rhoRoe * (dvel[i] - (*normal)[i] * dV) );
    dot_prod += velRoe[i] * dvel[i];
  }
  df1[4]   = std::abs(VnRoe) * ( 0.5 * qSqRoe * drhoTerm +
             rhoRoe * ( dot_prod - VnRoe * dV ) );
  // Calculate df2
  df2[0]   = 0.5 * std::abs(VnRoe + aRoe) *
             ( dp / ( aRoe * aRoe ) + rhoRoe * dV / aRoe );
  for( int i = 0 ; i < 3 ; ++i )
    df2[i+1] = df2[0] * ( velRoe[i] + (*normal)[i] * aRoe );
  df2[4] = df2[0] * ( h0Roe + VnRoe * aRoe );
  // Calculate df3
  df3[0]   = 0.5 * std::abs(VnRoe - aRoe) *
             ( dp / ( aRoe * aRoe ) - rhoRoe * dV / aRoe );
  for( int i = 0 ; i < 3 ; ++i )
    df3[i+1] = df3[0] * ( velRoe[i] - (*normal)[i] * aRoe );
  df3[4]   = df3[0] * ( h0Roe - VnRoe * aRoe );
  // Sum fluxes 0.5 * (fL + fR - df1 - df2 -df3 ) to yield flux at the interface
  (*mass) = 0.5 * ( fL[0] + fR[0] - df1[0] - df2[0] - df3[0] );
  for( int i = 0 ; i < 3 ; ++i )
    (*mom)[i] = 0.5 * ( fL[i+1] + fR[i+1] - df1[i+1] - df2[i+1] - df3[i+1] );
  (*energy) = 0.5 * ( fL[4] + fR[4] - df1[4] - df2[4] - df3[4] );

  return std::fabs( VnRoe ) + aRoe;
}

/// Van-leer flux vector splitting scheme implemeted by 
scalar van_leer(scalar *rhoL , vector *uL , scalar *pL ,
          scalar *rhoR , vector *uR , scalar *pR ,
          scalar *mass , vector *mom , scalar *energy ,
          vector *normal){

  scalar h0L, h0R,  VnL, VnR, uL2 , uR2,  aL, aR, macL, macR,eL,eR ,fL[5], fR[5] ,mdotL, mdotR;

  //  Left and Right Flow values
  uL2 = (*uL) & (*uL);
  uR2 = (*uR) & (*uR);
  h0L = gbygm1 * (*pL) / (*rhoL) + 0.50 * ( uL2 );
  h0R = gbygm1 * (*pR) / (*rhoR) + 0.50 * ( uR2 );
  VnL = (*uL) & (*normal);
  VnR = (*uR) & (*normal);
  aL  = std::sqrt(gama*(*pL)/(*rhoL));
  aR  = std::sqrt(gama*(*pR)/(*rhoR));
  macL = VnL/aL;
  macR = VnR/aR;
  eL   = (*pL)/(gm1) + 0.50*(*rhoL)*(uL2);
  eR   = (*pR)/(gm1) + 0.50*(*rhoR)*(uR2);
  mdotL = 0.25*((*rhoL)*aL*pow((1+macL),2));
  mdotR = -0.25*((*rhoR)*aR*pow((1-macR),2));
 
 ///  To compute the left face  fluxes
  if(macL<=-1.0){
    for (int i =0; i<5; i++){
      fL[i]= 0.0;
      }}
   if(macL>=1.0){
     fL[0] = (*rhoL) *VnL;
     fL[1] = (*rhoL)*(*uL)[0]*VnL+(*pL)*(*normal)[0];
     fL[2] = (*rhoL)*(*uL)[1]*VnL+(*pL)*(*normal)[1];
     fL[3] = (*rhoL)*(*uL)[2]*VnL+(*pL)*(*normal)[2];
     fL[4] = (eL+(*pL))*VnL;
   }
   if((macL>-1.0)&&(macL<1.0)){
     fL[0] = mdotL;
     fL[1] = mdotL*((*uL)[0]+(2-macL)*(aL)*(*normal)[0]/gama);
     fL[2] = mdotL*((*uL)[1]+(2-macL)*(aL)*(*normal)[1]/gama);
     fL[3] = mdotL*((*uL)[2]+(2-macL)*(aL)*(*normal)[2]/gama);
     fL[4] = mdotL*h0L;
     }

 ///  To compute the right face  fluxes
  if(macR>=1.0){
    for (int i =0; i<5; i++){
      fR[i]= 0.0;
      }}
   if(macR<=-1.0){
     fR[0] = (*rhoR) *VnR;
     fR[1] = (*rhoR)*(*uR)[0]*VnR+(*pR)*(*normal)[0];
     fR[2] = (*rhoR)*(*uR)[1]*VnR+(*pR)*(*normal)[1];
     fR[3] = (*rhoR)*(*uR)[2]*VnR+(*pR)*(*normal)[2];
     fR[4] = (eR+(*pR))*VnR;
   }
   if((macR>-1.0)&&(macR<1.0)){
     fR[0] = mdotR;
     fR[1] = mdotR*((*uR)[0]-(2+macR)*(aR)*(*normal)[0]/gama);
     fR[2] = mdotR*((*uR)[1]-(2+macR)*(aR)*(*normal)[1]/gama);
     fR[3] = mdotR*((*uR)[2]-(2+macR)*(aR)*(*normal)[2]/gama);
     fR[4] = mdotR*h0R;
   }

  //  fluxes  fL + fR to yield flux at the interface
  (*mass) =   fL[0] + fR[0] ;
  for( int i = 1 ; i < 4 ; ++i )
    (*mom)[i-1] =  fL[i] + fR[i] ;
  (*energy) =  fL[4] + fR[4] ;

  return 0.5 * ( std::fabs( VnL + VnR ) + aL + aR );
}   /// End of Vanleer scheme


