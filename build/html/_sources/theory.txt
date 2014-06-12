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

As the name suggests, these boundary conditions are imposed when the Mach number of the the flow entering or leaving the boundary of the computational domain is fully supersonic. Therefore, there is no influence of the downwind disturbances. One can make use of this property and prescribe the inflow/outflow boundary conditions. For a supersonic inflow boundary face, the flow parameters at the inlet are used to obtain the fluxes. For a supersonic outflow boundary face, the solution state extrapolated from within the interior computational domain is used for calculating the fluxes across that face.

  2. Slip-wall

The slip wall boundary condition physically imposes a zero mass flux crossing the rigid wall. This can be written mathematically as :math:`u_n = 0`. This immediately suggests the following numerical flux formulation for the wall boundary face,

.. math::
  F_n^{wall} = \left[
  \begin{array}{c}
  \rho u_n \\
  \rho u u_n + p n_x \\
  \rho v u_n + p n_y \\
  \rho w u_n + p n_z \\
  (e + p ) u_n
  \end{array} 
  \right]
  = \left[
  \begin{array}{c}
  0 \\
  p n_x \\
  p n_y \\
  p n_z \\
  0
  \end{array} 
  \right]
  :label: wall_bc

The pressure value at the wall boundary face is usually computed by extrapolation from the interior computational domain. The nearest cell centroid value extrapolation has been implemented in the present work.

  3. Riemann Extrapolation

The Riemann invariants for the 1D Euler's equation (isentropic flow assumption) can be utilized as far-field boundary conditions by making a 1D flow assumption normal to the boundary face. If the far-field is sufficiently away from wall boundaries isentropic flow assumption is valid. The invariants :math:`r_1,r_2` and :math:`r_3` normal to a given finite volume interface is defined in equations shown below,

.. math::
  r_1 = u_n + \frac{2 a }{\gamma - 1}, \hspace{5pt} r_2 = u_n + \frac{2 a }{\gamma - 1},
  \hspace{5pt} and \hspace{5pt} r_3 = \frac{p}{ \rho^{\gamma} }

One can detect the direction of propagation at the far-field boundary based on the sign of the Eigenvalues, which are defined as,

.. math::
  \lambda_1 = u_n + a, \hspace{5pt} \lambda_2 = u_n - a,
  \hspace{5pt} and \hspace{5pt} \lambda_3 = u_n

where, :math:`u_n` is the normal velocity to the interface and :math:`a` is the speed of sound at the interface. Based on the signs of the two variables four different cases can be considered as follows,

  * Supersonic Outflow (:math:`u_n \geq 0` and :math:`|u_n| \geq a`)

The boundary condition is exactly same as that in previous section.

  * Supersonic Inflow (:math:`u_n < 0` and :math:`|u_n| \geq a`)}

The boundary condition is exactly same as that in previous section.

  * Subsonic Outflow (:math:`u_n \geq 0` and :math:`|u_n| < a`)

In this case, :math:`\lambda_1 > 0, \hspace{5pt} \lambda_2 < 0 \hspace{5pt} and \hspace{5pt} \lambda_3 > 0`, which means that two characteristic waves are outgoing and one wave is incoming. The boundary values are then determined as follows,

.. math::
  :nowrap:

  \begin{eqnarray}
  u_n + a > 0 &\Longrightarrow& r_{1f} = r_{1e}\nonumber \\
  u_n - a < 0 &\Longrightarrow& r_{2f} = r_{2 \infty }\\
  u_n > 0 &\Longrightarrow& r_{3f} = r_{3e} \nonumber \\
  u_{nf} &=& \frac{ r_{1e} - r_{2 \infty} }{ 2 } \\
  a_f &=& \frac{\gamma - 1 }{4} ( r_{1e} - r_{2\infty}) \\
  u_{\parallel f,1} &=& u_{\parallel e,1} \\
  u_{\parallel f,2} &=& u_{\parallel e,2} \\
  \rho_f &=& \rho_e \\
  a_f &=& \sqrt{ \frac{\gamma p_f }{ \rho_f} }
  \end{eqnarray}

Note that the subscripts :math:`n, \hspace{5pt} \parallel, \hspace{5pt} f, \hspace{5pt} e` and :math:`\infty` denote values normal/parallel to the boundary face, boundary cell values (unknowns), extrapolated value from interior and free-stream value respectively.

  * Subsonic Inflow (:math:`u_n < 0` and :math:`|u_n| < a`)

In this case, :math:`\lambda_1 > 0, \hspace{5pt} \lambda_2 < 0` and :math:`\lambda_3 < 0`, which means that one characteristic wave is outgoing and two waves are incoming. The boundary values are then determined as follows,

.. math::
  :nowrap:

  \begin{eqnarray}
  u_n + a > 0 &\Longrightarrow& r_{1f} = r_{1e}\nonumber \\
  u_n - a < 0 &\Longrightarrow& r_{2f} = r_{2 \infty }\\
  u_n < 0 &\Longrightarrow& r_{3f} = r_{3 \infty} \nonumber \\
  u_{nf} &=& \frac{ r_{1e} - r_{2 \infty} }{ 2 } \\
  a_f &=& \frac{\gamma - 1 }{4} ( r_{1e} - r_{2\infty}) \\
  u_{\parallel f,1} &=& u_{\parallel e,1} \\
  u_{\parallel f,2} &=& u_{\parallel e,2} \\
  \rho_f &=& \rho_{\infty} \\
  a_f &=& \sqrt{ \frac{\gamma p_f }{ \rho_f} }
  \end{eqnarray}


4. Implicit non-reflecting boundary condition

Implicit boundary conditions are quite necessary for implementing Jacobian Free Newton-Krylov (JFNK) methods. Unlike the previous case, implicit non-reflecting boundary conditions fixes the normal boundary fluxes and not the boundary values. This way one ensures that the effect of boundary is implied in the residual calculation. Here one extrapolates the interior values for the left state and fix free-stream conditions as the right state and invokes the approximate Riemann flux solver to obtain the interfacial boundary fluxes and they are summed up to the neighboring cell residual. The implicit boundary described above is a weak boundary condition and hence is preferred more than the Riemann extrapolation due to the improved stability.


