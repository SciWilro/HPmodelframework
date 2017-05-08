# Prediction (regression + classification) of hybrid performance

Code and files for reproducing the analysis in de Abreu e Lima et al. (2017).

Note that the package 'caret' has many dependencies. Allow your machine for installation of additional packages when prompted. I highly recommend using R Open or other parallel-based computational methods, it might take 1-2 days to run the models.

script.R - contains all the code necessary to reproduce my results

metabolicData.txt - TSV input file worked on in 'script.R', contains all metabolic profiles (n = 1400) from hybrids and inbreds

biomassData.txt - TSV input file worked on in 'script.R', contains the biomass (FW) values per hybrid (n = 392) as well as the corresponding trial assignments ('2010' and '2012').

Please NOTE that the apparent difference in the sample sizes (i.e. 1400 and 392) is because these are raw data; in the first lines of the script you will be filtering profiles (i.e. hybrids in both biomass and metabolic data sets, profiles of crosses with insufficient replicates, etc).

For any inquiries:
francisco.lima278@gmail.com
:)
