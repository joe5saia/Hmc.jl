# Contains model estimation results

## official
This directory contains data for the model estimated with the official CPI and no noise in the data. 
The `*_summary.csv` files contain the means of each parameter averaged across the draws. 

## signals_official_noise_*_allsignal
This directory contains data for the model estimated with the official CPI and noise random, indpendent noise added to each observation in the data. 
The number in the directory name tells the relative std of the noise. 
The `*_summary.csv` files contain the means of each parameter averaged across the draws for each simulated sample. 
The `*_dispersion.csv` files contain the mean and std of each parameter averaged across the postior mean for each simulated sample