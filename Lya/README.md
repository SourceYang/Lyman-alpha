This is the main repository for Yuanyuan’s summer 2020 research project focusing on models of polarized Lyman alpha emission from reionization.

Before getting start, please use Chenxiao’s code (https://github.com/frankelzeng/reionization) $\textbf{gif\_190127\_2.c}$ to set up the background temperature by
```
cc -std=c99 -Wall -o gif gif_190127_2.c
```
or
```
gcc gif_190127_2.c -o gif -lm
```
and then
```
./gif T U > output_gif.txt
```
replacing T with blackbody temperature in K, U with front velocity in cm/s, and output_gif.txt with output file name

$\textbf{Lya.c}$ is a whole code including all parts of the model. It will generate a bunch of Lyman alpha photons and will propagate and scatter the photon until it escapes from the front, and then record its final status via
```
gcc Lya.c -o lya -lm
./lya T U y1H1 output_gif.txt n_photon seed > output_lya.txt
```
The corresponding inputs are the the source balckbody temperature T in K, front velocity U in cm/s, hydrogen density y1H1 in #/cm^-3, background file name (the output from gif_190127_2.c), number of photons n_photon, a random seed and the output file name

To test, you can use the default T = 5e4 K, the default U = 5e8 cm/s and the default y1H1 = 1.37e-4 cm^-3, if you want to change any of the parameter, please change the corresponding parameter when set up the background temperature as well
