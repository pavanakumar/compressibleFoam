   /// Register cmd line inputs
   Foam::argList::validOptions.set( "mach", "Mach number" );
   Foam::argList::validOptions.set( "aoa", "Angle of attack" );
   Foam::argList::validOptions.set( "cfl", "CFL number" );
 
   /// Default values
   scalar M_inf, aoa, CFL;

   if( args.optionReadIfPresent( "mach", M_inf ) == false )
     Info << "Need a Mach number input " << Foam::FatalError;
   /// read AOA
   if( args.optionReadIfPresent( "aoa", aoa ) == false ) aoa = 0.0;
   /// Read CFL
   if( args.optionReadIfPresent( "cfl", CFL ) == false ) CFL = 0.8;

