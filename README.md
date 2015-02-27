# PNASForageFish
Code and data used for PNAS paper "Fishing amplifies forage fish population collapses"

# Data
The main data file is allforagedata.csv. It contains all time series of total biomass (TB), spawning stock biomass (SSB), total number (
(TN), Recruitment (R), total catch (TC, total landings (TL), instantaneous fishing mortality rate (F), and exploition rate (ER).  
It also contains the stock ID codes, species, and long description for each stock.  

# Code
There are several R files, each prepared to run independent of the others.  While this makes for some long files, this 
allows the user to recreate analyses in our manuscript exactly and to explore variations.
ALL files use the data file allforagedata.csv.  All files therefore require you to specify the path to the folder 
where you hold this file
- Compare Collapse Frequences.R: looks at the year of collapse of stock, assigns to a decade, and then asks whether the frequency of 
collapse is different by decade, region
- Compare production and exploitation collapsed and non-collapsed stocks.R: Runs the analysis that compares mean surplus production and 
mean exploitation rate immediately before collapse or minimum biomass levels
- figure 2 analysis - runs the data analysis used to generate figure 2
- Figure 3 analysis - runs the randomization test (1,000 iterations) and plots the results.  This may take a long time to run
- Test Hypothetical Harvest Control Rule.R - performs the thought experiment described in manuscript in which fishing is suspended
whenever biomass is less than one half the long term mean

Other code include:
- Summarize B and Fdata.R: not needed, but is used to generate the .Rdata file that takes information in the "allforagedata.csv" file
and organizes it into matrices of biomass and exploitation rate.  This same code appears in several of the files above
