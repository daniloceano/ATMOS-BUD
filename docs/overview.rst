.. image:: _static/images/logo.jpg
   :alt: ATMOS-BUD logo
   :align: center


Overview
========

ATMOS-BUD is an advanced software suite designed for calculating heat, vorticity, and moisture balances within confined atmospheric regions. Originally developed for the Synoptic Meteorology III course at the Institute of Astronomy, Geophysics, and Atmospheric Sciences at the University of São Paulo, ATMOS-BUD serves as a key tool for both educational and research purposes. It is widely used by undergraduate and graduate students, as well as experienced researchers in atmospheric science.

Developed by doctoral student and PAE assistant Danilo Couto de Souza, under the mentorship of Professors Pedro Leite da Silva Dias and Ricardo Hallak, ATMOS-BUD embodies the fusion of academic rigor and practical application. The software efficiently processes atmospheric model and reanalysis data in NetCDF format, producing reliable outputs, including:
- CSV files containing spatial averages of each variable involved in the equations being analyzed, across different vertical levels and time intervals.
- A comprehensive NetCDF file that includes the spatial results for each calculated term, across all vertical levels and temporal snapshots.

ATMOS-BUD is designed with flexibility in mind, offering three distinct operational frameworks to accommodate various atmospheric phenomena:

1. **Fixed Domain:** Designed for studying slowly evolving atmospheric systems, such as convergence zones, where the spatial domain remains unchanged over time.
2. **Semi-Lagrangian Domain:** Ideal for dynamic analyses of moving systems like cyclones, this framework allows the spatial domain to follow the system's movement at each time step.
3. **Interactive Domain:** Providing the highest level of flexibility, this framework enables users to interactively define their domain of interest at any given time, making it a powerful tool for detailed, user-driven investigations.

ATMOS-BUD is not just a software tool; it is a gateway to understanding complex atmospheric dynamics, offering users an intuitive experience for both educational exploration and advanced research.

Thermodynamic Equation
-----------------------------------------

.. math::

    \frac{\partial T}{\partial t} = & -\mathbf{V}_h \cdot \nabla_h T & \text{(Horizontal Temperature Advection)} \\
                                    & + S_p \omega & \text{(Total Vertical Motion Effect)} \\
                                    & + Q & \text{(Diabatic Heating)}

where :math:`S_p` is an approximation of the static stability term, given by:

.. math::

    S_p = \frac{R T}{c_p p} - \frac{\partial T}{\partial p} = - \frac{T}{\theta} \frac{\partial \theta}{\partial p}

Vorticity Budget Equation
-------------------------------------------

.. math::

   \frac{\partial \zeta}{\partial t} = & -\mathbf{V}_h \cdot \nabla_h \zeta & \text{(Horizontal Advection)} \\
                                       & - \omega \frac{\partial \zeta}{\partial p} & \text{(Vertical Advection)} \\
                                       & - \beta v & \text{(Beta Term)} \\
                                       & - \zeta \nabla \cdot \mathbf{V}_h & \text{(Stretching Term)} \\
                                       & - f \nabla \cdot \mathbf{V}_h & \text{(Divergence Term)} \\
                                       & + \left( \frac{\partial \zeta}{\partial y} \frac{\partial u}{\partial p} - \frac{\partial \zeta}{\partial x} \frac{\partial v}{\partial p} \right) & \text{(Tilting Term)}

Water Budget Equation
---------------------

.. math::

   \int \frac{\partial q}{\partial t} \, dp = & -\int \mathbf{V}_h \cdot \nabla q \, dp & \text{(Integrated Moisture Advection)} \\
                                              & + \int S \, dp & \text{(Net Source/Sink Term)}
