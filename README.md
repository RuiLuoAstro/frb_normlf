# frb_normlf
FRB mock data simulator and a Bayesian MCMC framework to measure FRB luminosity function

## NOTES
If you would like to use this package to do FRB statistics, please cite the paper [Luo et al. 2018a](http://adsabs.harvard.edu/abs/2018MNRAS.481.2320L)

## Dependancies
**Python 2.7.x**: NOT 3.x.x, unless you would like to transform the cod in 2.7.x to 3.x.x 

**Numpy**: version 1.14.5 at least

**Scipy**: version 1.0.0 at least

**Matplotlib**

**PyMultiNest** (see https://github.com/JohannesBuchner/PyMultiNest for more details)

## Simulate mock FRB data
Examples: 

$ ./simufrb.py -alpha alpha -logls logls -ns Nfrb -thre flux_thre -dnu specwidth -type galaxy_type -o simfrb.txt

Options: 

-alpha  &emsp;&emsp;&emsp;&emsp;  Inputs the power-law index of FRB luminosity function

-logls  &emsp;&emsp;&emsp;&emsp;  Inputs the expoential cut-off of FRB luminosity function in logarithmic erg/s

-ns  &emsp;&emsp;&emsp;&emsp;  Inputs the FRB number you want to generate

-thre  &emsp;&emsp;&emsp;&emsp;  Sets Flux threshold in units of Jy

-dnu  &emsp;&emsp;&emsp;&emsp;  Sets the secptral width in units of MHz

-type &emsp;&emsp;&emsp;&emsp; Chooses galaxy case among ETG_NE2001, ETG_YMW16, LTG_NE2001, LTG_YMW16, ALG_NE2001, ALG_YMW16 (see [Luo et al. 2018a](http://adsabs.harvard.edu/abs/2018MNRAS.481.2320L) for more details)

-o &emsp;&emsp;&emsp;&emsp;  Outputs the mock data in .txt format

## Verify the mock data using PyMultiNest
*Use*:

$ ./run_simu.sh &emsp;&emsp;&emsp;&emsp;
Notes: better implement it in cluster where MPI was installed well

$ ./draw_sim.sh &emsp;&emsp;&emsp;&emsp;
Plot the posterior distribution contours of the mock data

## Measure the FRB LF of current sample
*Use*:

$ ./run_samp.sh &emsp;&emsp;&emsp;&emsp;
Notes: better implement it in the cluster where MPI was installed well

$ ./draw_samp.sh &emsp;&emsp;&emsp;&emsp;
Plot the posterior distribution contours of the real FRB sample

$ .condat.sh &emsp;&emsp;&emsp;&emsp;
Get the contour data, which contains the best inferred value and error area of each parameter. 
