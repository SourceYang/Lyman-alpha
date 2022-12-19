This is a follow-up project of model of polarized lyman-alpha emission from cosmic ionization front. This project takes information from the single front lyman alpha intensity calculation and applies it to cosmological simulations from 21cmFAST to calculate a realistic power spectrum of intensity and polarized intensity from reionization. 

To start, this assumes the user has already constructed a 21cmFAST cosmological simulation. Specifically, one needs an ionization box with xH_box and Gamma12_box information stored. 

EdgeFinding3d072622.py takes the 21cmFAST simulation and breaks the fronts into a collection of triangles. To accoplish this, the box is divided into 2D slabs and line segments ("pairs") are recorded that run along the half-ionized border of ionization bubbles. Next, "lines" are stored that connect "pairs between two adjacent slabs. Finally, "fronts" are calculated by building triangles from the lines and pairs, and information about the area, center, and direction of the front are recorded. From this, the location that Gamma12 should be calculated at for each triangle front is determined, and the Gamma12 values are recorded in a list. 

Also necessary for a power spectrum is an interpolation scheme ...

Then need to do calculation of power spectrum ...
