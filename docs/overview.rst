.. image:: _static/images/logo.jpg
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

Thermodynamic Equation
-----------------------------------------

.. math::

    \frac{\partial T}{\partial t} = & -\mathbf{V}_h \cdot \nabla_h T & \text{(Horizontal Temperature Advection)} \\
                                    & - \Sigma \omega & \text{(Vertical Temperature Advection)} \\
                                    & + Q & \text{(Diabatic Heating)}

Vorticity Budget Equation
-------------------------------------------

.. math::

   \frac{\partial \zeta}{\partial t} = & -\mathbf{V}_h \cdot \nabla_h \zeta & \text{(Horizontal Advection)} \\
                                       & - \frac{\zeta \omega}{\partial p} & \text{(Vertical Advection)} \\
                                       & - \beta v & \text{(Beta Term)} \\
                                       & - \zeta \nabla \cdot \mathbf{V}_h & \text{(Stretching Term)} \\
                                       & - f \nabla \cdot \mathbf{V}_h & \text{(Divergence Term)} \\
                                       & + \left( \frac{\partial \zeta}{\partial y} \frac{\partial u}{\partial p} - \frac{\partial \zeta}{\partial x} \frac{\partial v}{\partial p} \right) & \text{(Tilting Term)}

Water Budget Equation
---------------------

.. math::

   \int \frac{\partial q}{\partial t} \, dp = & -\int \mathbf{V}_h \cdot \nabla q \, dp & \text{(Integrated Moisture Advection)} \\
                                              & + \int S \, dp & \text{(Net Source/Sink Term)}
