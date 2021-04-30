# Lyman-alpha

This is the main repository for Yuanyuan’s summer 2020 research project focusing on models of polarized Lyman alpha emission from reionization.

Before getting start, please add r64c.h random number generator, and use Chenxiao’s code (https://github.com/frankelzeng/reionization) gif_190127_2.c to set up the background temperature, and please use outputs.txt as the output file (or will need to change the filename inside each code)

Lya.c is a whole code including all parts of the model. It will generate a bunch of Lyman alpha photons and will propagate and scatter the photon until it escapes from the front, and then record its final status. The inputs are the source balckbody temperature T, front velocity U, hydrogen density y1H1, number of photons n_photon and a random seed

The default T is 5e4 K, the default U is 5e8 cm/s and the default y1H1 is 1.37e-4 cm^-3, if you want to change any of the parameter, please change the corresponding parameter when set up the background temperature as well
