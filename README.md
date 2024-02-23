# AUC diabetes
This script helps calculate cummulated AUC at different timestamps.

The input is a dataframe with a column consisting of time-points and outcome-points. 
These pairs are used to calculate the AUC at various time points.
The output is a dataframe in long format.


The melt/pivot function converts the long format to a wide format.
