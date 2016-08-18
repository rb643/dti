# Still a work in progress...

So far contains 3 scripts

- [dti.m] (https://github.com/rb643/dti/blob/master/dti.m), this is the main script to run all functions. Just type help dti and it will tell you what to do and what your options are
- [runStats.m] (https://github.com/rb643/dti/blob/master/runStats.m), this runs some statistics on the output you would get from running dti.m. You can specify which output you want and which metric. Currently only runs a standard ANOVA + post-hoc permutation on global metrics.

Requirements:
- Matlab 2014B or higher
- Brain Connectivity Toolbox functions
- [getTop.m] (https://github.com/rb643/dti/blob/master/getTop.m), this is the script/function you need to extract the top set of nodes based on the degree

Note: so far only tested on Mac iOS El Capitan
