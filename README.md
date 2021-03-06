## Effects of landcover on avian community composition and co-abundance variation

This `BBS.occurrences` repository stores functions and datasets for a collaborative project aiming to analyse variation in avian community composition and interspecific 'interactions' by applying Conditional Random Fields (Clark et al. 2018) to count data using [North American Breeding Bird Survey](https://www2.usgs.gov/science/cite-view.php?cite=17) data.   
  
Code in the `Appendix_S1_CompileBBS` script (within the `Clark_etal_analysis` folder) was used for compiling BBS observations from years 2003 - 2009.   
  
`Appendix_S2_CompileLandcover` outlines how functions in the repository https://github.com/nicholasjclark/LandcoverMODIS were used for accessing 40 x 40km resolution landcover data from MODIS for each site in each year.  
  
`Appendix_S3_FlywayGrouping` demonstrates how we group observation sites into recognised migratory flyways and by counties using `R`'s ever-increasing GIS capabilities.  
  
`Appendix_S4_PhyloTraitData` describes how avian species' phylogenetic and functional
trait datasets were accessed to create pairwise distance matrices.  
  
`Appendix_S5_ClimateAltitudeData` describes how to download climate data from the [PRISM dataset](https://climatedataguide.ucar.edu/climate-data/prism-high-resolution-spatial-climate-data-united-states-maxmin-temp-dewpoint) and altitude data from [NASA](https://asterweb.jpl.nasa.gov/gdem.asp) for each sample observation. NDVI data was downloaded at 40km resolution using functions in the `gimms` package    
  
`Appendix_S6_StormData` demonstrates how to download storm data from the [NOAA storm events database](https://d396qusza40orc.cloudfront.net/repdata%2Fdata%2FStormData.csv.bz2).  
  
`Appendix_S7_SpeciesTraits` shows functions used to scrape BirdLife International for information on species' migratory status, create binary indicators of species' nesting behaviour and broad habitat category use using data in the [Land Condition Trend Analysis Avian Database](https://ecologicaldata.org/wiki/land-condition-trend-analysis-avian-database-ecological-guild-based-summaries), calculate Shannon Diversity indices of habitat and diet diversity for each species using data in the [EltonTraits database](https://figshare.com/articles/Data_Paper_Data_Paper/3559887) and extract information on average clutch size from the [Avian Life History Dataset](http://www.esapubs.org/archive/ecol/E088/096/default.htm).
  
## References
Clark, N.J., Wells, K., Lindberg, O. (2018). Unravelling changing interspecific interactions across environmental gradients using Markov random fields. [*Ecology*](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2221) DOI: https://doi.org/10.1002/ecy.2221&nbsp;|&nbsp;[PDF](http://nicholasjclark.weebly.com/uploads/4/4/9/4/44946407/clark_et_al-2018-ecology.pdf)  
  
Poisot, T., Canard, E., Mouillot, D., Mouquet, N. & Gravel, D. (2012). The dissimilarity of species interaction networks. [*Ecology Letters*](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12002), 15, 1353-1361.


*This project is licensed under the terms of the MIT + file LICENSE*

