# bilinear isotropic hardening for 2D plane stress 

This code provides a 2D MATLAB example demonstrating bilinear (linear) isotropic hardening under plane stress conditions at a single Gauss point. Several options for strain excitation are included:
1. `normal_x`: strain excitation applied in the x-direction ($\varepsilon_\mathrm{xx} \neq 0$, $\varepsilon_\mathrm{yy} = 0$, $\gamma_\mathrm{xy} = 0$)
2. `normal_y`: strain excitation applied in the y-direction ($\varepsilon_\mathrm{xx} = 0$, $\varepsilon_\mathrm{yy} \neq 0$, $\gamma_\mathrm{xy} = 0$)
3. `shear`: strain excitation applied in the xy-direction ($\varepsilon_\mathrm{xx} = 0$, $\varepsilon_\mathrm{yy} = 0$, $\gamma_\mathrm{xy} \neq 0$)
4. `combined`: combined strain excitation applied in the x, y, and xy-direction ($\varepsilon_\mathrm{xx} \neq 0$, $\varepsilon_\mathrm{yy} \neq 0$, $\gamma_\mathrm{xy} \neq 0$)

The C++ version of this code has been implemented in GetDP (ONELAB) to simulate the elasto-plastic behavior of electrical machines. Nevertheless, it is also possible to calculate other geometries or models. 

The following workflow describes the algorithm functionality:
[![workflow]](images/QR_Codesvg.svg)
