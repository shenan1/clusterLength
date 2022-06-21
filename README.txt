turing-wind
===========

Turing energy project (wind downscaling to sub-hourly values)

This README file is for R code, clusterLength.R, written by Shenan Grossberg s.grossberg@exeter.ac.uk

Purpose of the code
===================

The purpose of the code is to calculate cluster lengths of wind ramps and wind droughts for two modelled and two observed data sets. The cluster lengths are calculated for the whole time period and are also categorised by season.

The results are visualised as boxplots, both with and without outliers, the latter to provide higher resolution of the inter-quartile range. The boxplots of the modelled data and observed data results are placed side by side to facilitate comparison.

How to run the code
===================

The code reads data files from and outputs plots to a directory, called sgrossberg_clusterLength, which is saved on the jceei-turing-wind S3 bucket. So that the paths work, download this directory and to save this code to it.

The code was written in R version 4.0.3. It was set up to run in RStudio, but could also be run from the command line. The libraries to install are listed in the first section of the code.

The first seven sections of the code need only be run once.

Section 8 prompts the user to enter the modelled and observed data sets that they want to compare. In practice, both modelled data sets should be compared with both observed data sets, so the code should be run 4 times to generate 4 comparisons.

The modelled and observed data frames are saved at the end of sections 9 and 10, respectively. So, if they have already been saved, they can be read in sections 11 and 12, respectively, instead of being recreated. Note that the code in section 8 should be run beforehand, to ensure the intended data sets are compared.

Having created (or read) the data frames, sections 13-21 should be run in numerical order.

There was some concern about the accuracy of the K13 data from 2019 onwards, so the results were also obtained only until the end of 2018. However, it wasn’t decided whether this data was inaccurate, so the lines in sections 12 and 17, which remove the data after 2018, have been commented out.

Note that section 22 is only for the EURO-CORDEX modelled data at the Europlatform site. It was written to investigate the extreme seasonal wind ramp length. It generates a time series plot of the wind ramp events in Spring 2004 and calculates the run length for this time period.
