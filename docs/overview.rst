.. image:: docs/_static/images/logo.jpg
   :alt: ATMOS-BUD logo
   :align: center


Overview
========

The ATMOS-BUD is a comprehensive software suite designed for the calculation of heat, vorticity, and moisture balances within limited areas of the atmosphere. Initially developed for the Synoptic Meteorology III course at the Institute of Astronomy, Geophysics, and Atmospheric Sciences of the University of São Paulo, ATMOS-BUD serves as a pivotal educational and research tool for undergraduate and graduate students, as well as seasoned researchers in the field of atmospheric sciences.

Crafted by doctoral student and PAE assistant Danilo Couto de Souza, under the expert guidance of professors Pedro Leite da Silva Dias and Ricardo Hallak, ATMOS-BUD stands as a testament to the intersection of academic rigor and practical application. The software adeptly processes atmospheric model and reanalysis data in NetCDF format, delivering robust outputs that include:

- CSV files encompassing spatial averages of each variable integral to the equations under study, across various vertical levels and time points.
- A comprehensive NetCDF file containing the spatial results of each calculated term, for every vertical level and temporal snapshot.

ATMOS-BUD is engineered with versatility in mind, offering three distinct operational frameworks to cater to a broad spectrum of atmospheric phenomena:

1. **Fixed Domain:** Tailored for the study of slowly evolving atmospheric systems, such as convergence zones, where spatial domains remain constant over time.
2. **Semi-Lagrangian Domain:** Ideal for the dynamic analysis of mobile systems like cyclones, this framework allows the spatial domain to move with the system of interest at each time step.
3. **Interactive Domain:** Providing the ultimate flexibility, this framework enables users to interactively select their domain of interest at any given time step, making it a powerful feature for detailed, user-guided investigations.

ATMOS-BUD is more than a software program—it's a gateway to understanding the complex dance of atmospheric dynamics, offering users a rich, intuitive experience for both educational exploration and advanced research.

Quasi-Geostrophic Thermodynamic Equation
-----------------------------------------

.. math::

    \begin{align*}
    \frac{\partial T}{\partial t} = & \underbrace{-\mathbf{V}_h \cdot \nabla_h T}_{\text{Horizontal Temperature Advection}} \\
    & \underbrace{- \Sigma \omega}_{\text{Vertical Temperature Advection}} \\
    & \underbrace{+ Q}_{\text{Diabatic Heating}}
    \end{align*}


Quasi-Geostrophic Vorticity Budget Equation
-------------------------------------------

.. math::

   \frac{\partial \zeta}{\partial t} = & -\mathbf{V}_h \cdot \nabla_h \zeta \quad \text{(Horizontal Advection)} \\
   & - \frac{\zeta \omega}{\Delta p} \quad \text{(Vertical Advection)} \\
   & - \beta v \quad \text{(Beta Term)} \\
   & - \zeta \nabla \cdot \mathbf{V}_h \quad \text{(Stretching Term)} \\
   & - f \nabla \cdot \mathbf{V}_h \quad \text{(Divergence Term)} \\
   & + \left( \frac{\partial \Omega}{\partial y} \frac{\partial u}{\partial p} - \frac{\partial \Omega}{\partial x} \frac{\partial v}{\partial p} \right) \quad \text{(Tilting Term)}

Water Budget Equation
---------------------

.. math::

   \int \frac{\partial q}{\partial t} \, dp = & -\int \mathbf{V}_h \cdot \nabla q \, dp \quad \text{(Integrated Moisture Advection)} \\
   & + \int S \, dp \quad \text{(Net Source/Sink Term)}
