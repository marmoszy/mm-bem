## Jupyter notebooks with python mm-bem codes using bempp library


The notebooks show the folowing cases of acoustic scattering of plane wave in water:

1. soft scattering - [soft.ipynb](soft.ipynb) (without bempp library), [soft\_bempp.ipynb](soft_bempp.ipynb)
2. fluid scattering - [fluid\_bempp.ipynb](fluid_bempp.ipynb)
3. fluid with one inclusion scattering - [shell\_bempp.ipynb](shell_bempp.ipynb)
4. fluid with two inclusions scattering - [scat3\_bempp.ipynb](scat3_bempp.ipynb)

Each example contains:

- BEM python code
- mesh 3D preview
- scattering function polar plot in dB
- two-column scattering function (degree, scattering pressure)

Moreover it contains calculations of yellowfin tuna swimbladder from dorsal aspect using several technologies:
- scipy lu_solver - [yft\_scipy.ipynb](yft_scipy.ipynb)
- bempp gmres solver - [yft\_bempp.ipynb](yft_bempp.ipynb)
- pytorch solver - [yft\_pytorch.ipynb](yft_pytorch.ipynb)
- taichi + scipy lu_solver - [yft\_taichi.ipynb](yft_taichi.ipynb)

Examples:

1. Scattering from soft target
![](soft_bempp_1.svg)
![](soft_bempp_2.svg)

2. Scattering from fluid target
![](fluid_bempp_1.svg)
![](fluid_bempp_2.svg)

2. Scattering from fluid with inclusion
![](shell_bempp_1.svg)
![](shell_bempp_2.svg)

2. Scattering from fluid with two inclusions
![](scat3_bempp_1_0.svg)
![](scat3_bempp_2_0.svg)

Note:
All mesh datafiles are in [msh](../msh/) directory

MM 31.10.2026
