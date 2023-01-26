# cyclone_thermodynamics

This program solves the Quasi-Geostrophic Thermodynamic equation for a closed region on the atmosphere, explicitly computing each term and estimating the diabatic heating term (R_[t]) as a residual:

![image](https://user-images.githubusercontent.com/56005607/214878079-a359b897-2388-4197-bd95-3a0d0038ceda.png)

Where the first term in the right hand size is the temperature tendency, the second one is the horiztonal temperature advection and the last one is the vertcial velocity, in pressure levels, times the static stability term.
