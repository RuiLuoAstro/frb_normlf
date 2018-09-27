# frb_normlf
Measuring the FRB normalized luminosity function in a Bayesian MCMC framework

# Dependancies
Python 2.7

Numpy (version 1.14.5 at least)

Scipy (version 1.0.0 at least)

Matplotlib (version 2.2.2 at least)

PyMultiNest (see https://github.com/JohannesBuchner/PyMultiNest for more details)

# Simulate mock FRB data
examples: ./simufrb.py -alpha alpha -logls logls -ns Nfrb -thre flux_thre -dnu specwidth -type galaxy_type -o simfrb.txt

Options: 
-alpha          Input the power-law index of FRB luminosity function
