/// Perform non-blocking send/recv of all processor boundary data
/// duing the nth iteration so that the values are available to do the 
/// n+1 th time step
forAll ( mesh.boundaryMesh() , ipatch ) {
  word BCType = mesh.boundaryMesh().types()[ipatch];
  if( BCType == "processor" ) {
    rho.boundaryField()[ipatch].initEvaluate();
    U.boundaryField()[ipatch].initEvaluate();
    p.boundaryField()[ipatch].initEvaluate();
  } 
} 

