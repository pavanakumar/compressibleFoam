.. figure:: figures/logo_small.png 

Installation
=============

Cloning compressibleFoam source from github.com
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Make sure you have git version control system installed in your machine. Then
clone the compressibleFoam `repository <https://github.com/pavanakumar/compressibleFoam>`_ from github.com. 

  >>> git clone https://github.com/pavanakumar/compressibleFoam
  Cloning into 'compressibleFoam'...
  remote: Counting objects: 161, done.
  remote: Compressing objects: 100% (146/146), done.
  remote: Total 161 (delta 70), reused 54 (delta 12)
  Receiving objects: 100% (161/161), 957.67 KiB | 295.00 KiB/s, done.
  Resolving deltas: 100% (70/70), done.
  Checking connectivity... done.
  >>> cd compressibleFoam

Now checkout the master branch, which is the working branch.

  >>> git checkout master
 
Compiling compressibleFoam
^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the wmake utility from OpenFOAM to compile the compressibleFoam solver. Remember to source the OpenFOAM bashrc or cshrc script (which ever is applicable) before invoking wmake.
  
  >>> source $OF_INSTALL/etc/bashrc
  >>> wmake

