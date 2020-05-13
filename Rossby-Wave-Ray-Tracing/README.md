# Rossby-Wave-Ray-Tracing

This repo documents my scripts to conduct Rossby wave ray tracing under the framework of classic barotropical Rossby wave theory. All my efforts are based on 2 papers: Karoly, 1983 and Shaman 2005. These two papers can be found in the docs directory.

### Equations

$$
u_g = \bar{u}_M + [(k^2-l^2)\bar{q}_y-2kl\bar{q}_x]/K^4
$$

$$
v_g = \bar{v}_M + [2kl\bar{q}_y+(k^2-l^2)\bar{q}_x]/K^4
$$

$$
\dfrac{d_gk}{dt}=-k\bar{u}_x-l\bar{v}_x+(\bar{q}_{xy}k-\bar{q}_{xx}l)/K^4
$$

$$
\dfrac{d_gl}{dt}=-k\bar{u}_y-l\bar{v}_y+(\bar{q}_{yy}k-\bar{q}_{xy}l)/K^4
$$

$$
K^4=(k^2+l^2)^2
$$

$$
\omega=\bar{u}_Mk+\bar{v}_Ml+(\bar{q}_xl-\bar{q}_yk)/K^2
$$



In  these equations, several variables are determined by the background flow, which do not change over time but varies in space. These variables are: $\bar{u}_M, \bar{v}_M,q$ and their spatial derivatives.

### Instruction

Under `data`  directory, there is an simple NCL script to create input file.

Note that 90S and 90N should be excluded when create input file because 90S and 90N can not be dressed under Mercator projection.

`namelist.py`: some parameters

`util_func.py`: some utility functions

`trace.py`: this is the main part of ray tracing. Run it: `python trace.py`

After executing `trace.py`, the ray path file would be created under `output` directory.

Then you can modify `plot.ncl` accordingly and plot the ray path...





