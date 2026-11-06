## Jupyter notebooks with python mm-bem codes


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

Moreover it contains calculations of the target strength angular pattern from yellowfin tuna swimbladder model using several technologies:
- scipy lu_solver - [yft\_scipy.ipynb](yft_scipy.ipynb)
- bempp gmres solver - [yft\_bempp.ipynb](yft_bempp.ipynb)
- pytorch solver - [yft\_pytorch.ipynb](yft_pytorch.ipynb)
- taichi + scipy lu_solver - [yft\_taichi.ipynb](yft_taichi.ipynb)

![](out/yft_taichi_1.svg)

Fig. 1. Yellowfin tuna swimbladder model 

![](out/yft_taichi_3.svg)

Fig. 2. Yellowfin tuna swimbladder target strength angular pattern at 38kHz

<!--
*Tab. 1. The calculation times of scattering pattern for yellowfin tuna swimbladder model (fname = "YFT_swimbladder_origin.msh" ) containing 7502 vertices nad 15000 triangles (using vertices DOF=7502) on MacBookPro M1 with Sequoia 15.7.7.*

|<div align="center">YFT_swimbladder_origin.msh<br/>7502x15000</div>| DOF=7502 |
|:---|:---:|
|1. scipy|2m 41.5s|
|2. bempp|45.5s|
|3. pytorch|28.9s|
|4. taichi|7.2s|
-->

*Tab. 1. The calculation times of the target strength angular pattern for a yellowfin tuna reduced model (fname = "simplify10_YFT_swimbladder_origin.msh") containing 752 vertices nad 1500 triangles (using triangles centers DOF=1500) on MacBookPro M1 with Sequoia 15.7.7.*

|<div align="center">simplify10\_YFT\_swimbladder_origin.msh<br/>752x1500</div>|DOF=1500|
|:---|:---:|
|scipy|11.3s|
|bempp (gmres)|3m 20.8s|
|pytorch|20.0s|
|taichi|3.2s|

Scattering function examples:

1. Scattering from soft target
   
![](out/soft_bempp_1.svg)![](out/soft_bempp_2.svg)

2. Scattering from fluid target

![](out/fluid_bempp_1.svg)![](out/fluid_bempp_2.svg)

3. Scattering from fluid with inclusion

![](out/shell_bempp_1.svg)![](out/shell_bempp_2.svg)

4. Scattering from fluid with two inclusions

![](out/scat3_bempp_1_0.svg)![](out/scat3_bempp_2_0.svg)

Note:
All mesh datafiles are in [msh](../msh/) directory (datapath = "../msh/" variable)

MM 31.10.2026
