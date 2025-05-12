The file above is divided into 2 sections. The first section computes the Shannon Diversity Index for the years listed in variable years.keep.  
The second part of the code takes canopy cover percentage for each plot that has it, then calculates percent change over time by searching for the previous measurement (PREV_PLT_CN) and using the formula change = ((new - old)/|old|) * 100
