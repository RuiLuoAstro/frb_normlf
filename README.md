# frb_normlf

FRB mock data simulator and a Bayesian MCMC framework to measure FRB luminosity function

## NOTICE

If you would like to use this code to do FRB statistical studies, please cite the paper [Luo et al. 2018, MNRAS, 481, 2320](http://adsabs.harvard.edu/abs/2018MNRAS.481.2320L)

## Dependencies

### Python (2.7.x), Numpy (1.14 at least), Scipy (1.0.0 at least), PyMultiNest (see https://github.com/JohannesBuchner/PyMultiNest for more details), Matplotlib

## Simulate mock FRB data
### Example: 

$ ./simufrb.py -alpha alpha -logls logls -ns Nfrb -thre flux_thre -dnu specwidth -type galaxy_type -o simfrb.txt

### Options: 

-alpha  &emsp;&emsp;&emsp;&emsp;  **Inputs the power-law index of FRB luminosity function**

-logls  &emsp;&emsp;&emsp;&emsp;  **Inputs the expoential cut-off of FRB luminosity function in logarithmic erg/s**

-ns  &emsp;&emsp;&emsp;&emsp;&emsp;  **Inputs the FRB number you want to generate**

-thre  &emsp;&emsp;&emsp;&emsp;  **Sets Flux threshold in units of Jy**

-dnu  &emsp;&emsp;&emsp;&emsp;&emsp;  **Sets the secptral width in units of MHz**

-type &emsp;&emsp;&emsp;&emsp; **Chooses host galaxy case among ETG_NE2001, ETG_YMW16, LTG_NE2001, LTG_YMW16, ALG_NE2001, ALG_YMW16 (see [Luo et al. 2018a](http://adsabs.harvard.edu/abs/2018MNRAS.481.2320L) for more details)**

-o &emsp;&emsp;&emsp;&emsp;&emsp;  **Outputs the mock data in .txt format**

## Verify the mock data using PyMultiNest
### Use:

$ ./run_simu.sh &emsp;&emsp;&emsp;&emsp;
**Notes: better implement it in cluster where MPI was installed well**

$ ./draw_sim.sh &emsp;&emsp;&emsp;&emsp;
**Plot the posterior distribution contours of the mock data**

## Measure the FRB LF of current sample
### Use:

$ ./run_samp.sh &emsp;&emsp;&emsp;&emsp;
**Notes: better implement it in the cluster where MPI was installed well**

$ ./draw_samp.sh &emsp;&emsp;&emsp;&emsp;
**Plot the posterior distribution contours of the real FRB sample**

$ .condat.sh &emsp;&emsp;&emsp;&emsp;
**Get the contour data, which contains the best inferred value and error area of each parameter** 

## Appendix
### The use of mcmc_simu.py and mcmc_samp.py

**Examples**:

$ ./mcmc_simu.py -f inputfile -o outputfile -g galaxy_type

$ ./mcmc_frb.py -upper -f inputfile -o outputfile -g galaxy_type

**Options**:

-f &emsp;&emsp;&emsp;&emsp;&emsp; **Inputs the FRB catalog**

-o &emsp;&emsp;&emsp;&emsp;&emsp; **Outputs the posterior results with name of**

-g &emsp;&emsp;&emsp;&emsp;&emsp; **Choose host galaxy cases among ETG_NE2001, ETG_YMW16, LTG_NE2001, LTG_YMW16, ALG_NE2001, ALG_YMW16**

-upper &emsp;&emsp;&emsp; **Bool option, choose uniform prior for L0 or not, the other option is uniform prior for logL0**
