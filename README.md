# Lyman-alpha

This is the main repository for Yuanyuan’s summer 2020 research project focusing on models of polarized Lyman alpha emission from reionization.

Before getting start, please add ran2() random number generator, and use Chenxiao’s code (https://github.com/frankelzeng/reionization) gif_190127_2.c to set up the background temperature, and please use outputs.txt as the output file (or will need to change the filename inside each code)

There are several separate parts of the model:

test.c is the original version of the whole code

Lya.c is a whole code including all parts of the model. It will generate a bunch of Lyman alpha photons and will propagate and scatter the photon until it escapes from the front, and then record its final status. The inputs are the source balckbody temperature T, front velocity U, hydrogen density y1H1, number of photons n_photon and a random seed

The default T is 5e4 K, the default U is 5e8 cm/s and the default y1H1 is 1.37e-4 cm^-3

Other codes test each part of the main code separately:

Lya_position.c will generate a bunch of initial position for the Lyman alpha photons, and the inputs are T, U, y1H1, n_photon and a random seed

Lya_frequency.c will generate a bunch of initial frequency for the Lyman alpha photons, and the inputs are background temperature of hydrogen THI, n_photon and a random seed

Lya_propagate.c will calculate the propagating distance for pre-set Lyman alpha photons, and the inputs are T, U, y1H1, THI, frequency of the Lyman alpha photon L_v, n_photon and a random seed

Lya_scatter.c will calculate the status of pre-set Lyman alpha photons after emitting by hydrogen, and the inputs are THI, direction of the Lyman alpha photon L_mu, polarization of the Lyman alpha photon L_p, L_v, n_photon and a random seed

Finally, Lya.ipynb will draw the plot of the distribution for each output file
