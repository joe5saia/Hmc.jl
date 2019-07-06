#HMC.jl

This repo contains code to estimate the implementation of the Hiden Markov Chain model.

## File Structure

* `code` contains scripts to run stuff. Some of these are softlinked into the project root directory to
make it easier to run on the cluster

* `src` contains the code for the module `Hmc.jl` which contains all the functions used in  the estimation

* `test` contains unit test code

* `data` contains the data for the project
    -- `raw` contains the raw inflation data
    -- `output` contains output files from the estimation scripts

* `thoughts` contains latex code outline various procedures



## Data

The estimation code outputs the results from each Gibbs draw for the various parameters in the model.
If 1000 samples are drawn in the estimation for each date then the csv files will have 1000 lines for the
samples drawn for date range estimated. Each line starts with the end date of the estimation. For each
of these csv's there is also a second csv with the mean of the samples for each date. In this csv there is
exactly one line for each date. These csv's will have `_summary` as suffix.
The directory structure should make clear what's being estimated

* filtered_forecasts_x.csv are the 1-12 month ahead forecasts for inflation for the estimation ending on the given date
* filtered_means_x.csv are the means for the three states for the estimation ending on the given date
* filtered_variances_x.csv are the variances for the three states for the estimation ending on the given date
* filtered_state_probs_x.csv are the state probabilities for the three states for the last date of the estimation ending on the given date
* filtered_trans_probs_x.csv are the transition probabilites for the three states for the estimation ending on the given date. These are outputed by `A[:]` to transform the matrix to a vector and need to be reshaped in column major order (default in Julia and Matlab but not Python)
* smoothed_state_probs_x.csv are the state probabilities for the three states for each date of the sample from the last estimation



## Code
`run_hmm.jl` is the main script which estimates the HMM model. It takes a command line argument 
to determine which inflation series to use in the estimation. Call the script as either 
`julia --project=. -p <number of cores to use> run_hmm.jl official` or 
`julia --project=. -p <number of cores to use> run_hmm.jl alter`

`aggregate.jl` is a secondary script that calculats aggreate statistics

The root directory has shell scripts that can be submited to Slurm to run the code
