.. figure:: figures/logo_small.png 

compressibleFoam Theory
^^^^^^^^^^^^^^^^^^^^^^^^

Finite Volume Methods for Euler's Equation
-------------------------------------------

In the present chapter, we briefly derive the discretized finite volume scheme for the invscid compressible Euler's equations. The formulation is cell-centered with first order accuracy. The implementation of the same using OpenFOAM libraries is discussed. Finally, the parallelization of the solver using OpenFOAM library is discussed.

Governing Equations
--------------------
Shown in equations :eq:`comp_euler_cont`-:eq:`comp_euler_energy` are the compressible Euler equations in conservative form. The momentum equation :eq:`comp_euler_mom` is arranged such that the flow variables are on the left hand side and the driving potential terms are on the right. 

.. math::
  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
  :label: comp_euler_cont
.. math::
  \frac{\partial (\rho \mathbf{u})}{\partial t}+ \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \rho \mathbf{f}
  :label: comp_euler_mom
.. math::
  \frac{\partial (\rho e)}{\partial t}+ \nabla \cdot (\rho \mathbf{u} e ) = -\nabla \cdot (p \mathbf{u})
  :label: comp_euler_energy

The governing equations used for the CFD calculations are the compressible Euler's equations in Cartesian coordinates. The equations were introduced in differential form in equations :eq:`comp_euler_cont`- :eq:`comp_euler_energy`, but for the finite volume discretization, the integral form of the equations is employed. The integral equations in vector form is shown in equation :eq:`euler_src_integral`.

.. math::
  \int_V{\dot{Q} dV} + \int_S{\mathbf{n}\cdot \mathbf{F} }  = 0
  :label: euler_src_integral

where,

.. math::
  Q = \left[ \begin{array}{c}
  \rho \\
  \rho u \\
  \rho v \\
  \rho w \\
  e  \end{array} \right],
  \hspace{15pt}
  F_n = \mathbf{n \cdot F} = \left[ \begin{array}{c}
  \rho \bar{u} \\
  \rho u \bar{u} +p n_x\\
  \rho v \bar{u} +p n_y \\
  \rho w \bar{u} +p n_z \\
  e \bar{u}  + p u_n \end{array} \right]
  :label: dummy_1

The vector :math:`\mathbf{n}(n_x,n_y,n_z)` is the unit outward pointing normal of the surface of integration :math:`S`. The dot on top of any variable denotes differentiation with respect to time. The fluid velocity vector is :math:`\mathbf{u}(u,v,w)` and the grid velocity is denoted by :math:`\mathbf{v}`. The quantity :math:`u_n=\mathbf{n \cdot u}` and :math:`v_n=\mathbf{n \cdot v}` are the velocity components in the normal direction. Similarly, we have :math:`\bar{u}=\mathbf{n \cdot (u - v)}=u_n-v_n`.

Finite Volume Discretization
----------------------------

The volume integrals in equation :eq:`euler_src_integral` are evaluated using the cell center value :math:`Q_c` (piecewise constant) and the cell volume :math:`V_c` as shown in equation :eq:`approx_vol_int`.

.. math::
  \int_V{QdV} \approx Q_c V_c
  :label: approx_vol_int

.. _piecewise_const_reconst:

.. figure:: figures/piecewise_const_reconst.png

  Piecewise-constant reconstruction

In this piecewise-constant reconstruction one has to solve a local Riemann problem at every edge indicated by the fractional indices :math:`i+1/2$, $i-1/2` shown in figure :ref:`piecewise_const_reconst`. The face average flux :math:`F_n` (piecewise-constant) is approximated as the face center value and it is obtained by solving the Riemann problem at the cell edge (given the left and right states). Introducing the above approximations, the discretized equation can be written as shown in equation :eq:euler_src_integral_discrete.

.. math::
  \dot{Q_c} + \sum_{i}^{n}{\frac{F_n S_i}{V_i}}  = 0
  :label: euler_src_integral_discrete

The flux :math:`F_n` can be evaluated using an Approximate Riemann solver for example the Roe scheme, using the left and right states of a face :math:`Q_L` and :math:`Q_R` as shown in equation (:eq:roe_scheme). 

.. math::
  \hat{F}_{approx}(Q_L,Q_R,S)=\frac{1}{2}\left[\hat{F}(Q_L,S)+\hat{F}(Q_R,S) - \|\tilde{A}(Q_L,Q_R,S)\|(Q_R-Q_L)\right]
 :label: roe_scheme

The discrete equation shown in equation :eq:`euler_src_integral_discrete` can be rearranged and written in a semi-discrete form as an Ordinary Differential Equation (ODE) shown in equation :eq:`rk_lu`. The operator :math:`L` consists of the discrete flux operator and the contribution due to source term. The explicit Euler time integration method has been implemented to solve this ODE as shown in equation :eq:`euler_t`.

.. math::
  \dot{Q} = L[Q]
  :label: rk_lu

.. math::
  Q^{n+1} = Q^n + \frac{\Delta t}{v_{cell}} L[Q^n]
 :label: euler_t

where, :math:`\Delta t` is the time step and :math:`v_{cell}` is the cell volume. :math:`L[Q^n]` is the right hand side of the discrete equation (residue).

Local Time Stepping (LTS)
--------------------------

The solution is accelerated to steady state using the local time stepping (LTS). In equation :eq:`euler_t` the :math:`\frac{\Delta t}{v_{cell}}` term is replaced by the CFL number (:math:`N_c`) and local maximum eigenvalue integral as shown in equation :eq:`lts`.

.. math::
  \frac{\Delta t}{v_{cell}} = \frac{N_c}{ \displaystyle\sum_{i=1}^{N_f} \left( |U_n| + a \right) dS }
  :label: lts

where, :math:`N_f` is the total number of faces in a cell and :math:`dS` is the face area of the :math:`i^{th}` face.

Boundary Condition
-------------------

The following types of boundary conditions have been implemented, namely, 

  1. Supersonic Inflow
  2. Extrapolated Outflow
  3. Slip-wall
  4. Riemann Extrapolation
  5. Implicit non-reflecting boundary condition

