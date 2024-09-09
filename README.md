# transcontinental_common_garden

Readme file associated with the paper Albano et al., submitted to Ecological Monographs in September 2024. Title: Adaptation to climate in the native and introduced ranges of a cosmopolitan plant

In this paper, we examine the strength and environmental drivers of local adaptation in Trifolium repens in multiple ranges, along with the ability of T. repens to adapt in the face of environmental change. We use a large-scale trans-continental common garden experiment using 1858 individuals from 96 populations of T. repens transplanted into four experimental gardens, located in the north (Mississauga, ON, Canada) and south (Lafayette, LA, USA) of the introduced North American range, and the north (Uppsala, Sweden) and south (Montpellier, France) of the native European range. This experiment addressed five research questions: (Q1) Are T. repens populations locally adapted to their climate of origin and does the strength of local adaptation differ across space? (Q2) Does T. repens exhibit stronger local adaptation in its native range than in its introduced range? (Q3) What bioclimatic variables have the strongest impact on fitness at different latitudes in the native and introduced ranges of T. repens? (Q4) Are there adaptation lags to climate change in T. repens, such that populations originating from warmer and/or drier climates exhibit greater fitness than local populations? (Q5) What are the roles of cyanogenesis and herbivory in driving local adaptation? Three datafiles are used in order to address these research questions, with each of these files described below.

1. In order to assess local adaptation, we use a measure of climatic distance that was calculated using a Principal Component Analysis, combining six bioclimatic variables that were extracted from WorldClim (Version 2.1, February 2023) at a resolution of 30 seconds for each original T. repens population and common garden site (O’Donnell and Ignizio, 2012). The variables included were mean annual temperature (“BIO1”; MAT), maximum temperature of the warmest month (“BIO5”), minimum temperature of the coldest month (“BIO6”), annual precipitation (“BIO12”), precipitation of the wettest month (“BIO13”), and precipitation of the driest month (“BIO14”). Each variable was z-score standardized. These z-score standardized values for each population of origin and each common garden location can be found in Albano_et_al_WorldClim_z-standardized.csv

Columns A through F contain bioclimatic variables from BIO1 through BIO14 in increasing numerical order. Rows are not labeled because the code for conducting the PCA will not run properly if the dataframe does not contain only numerical values and column headers.
- Rows 2 through 48 contain bioclimatic data for all North American T. repens populations of origin
- Rows 49 through 97 contain bioclimatic data for all European T. repens populations of origin
- Row 98 contains bioclimatic data for the Lafayette, LA, USA common garden location
- Row 99 contains bioclimatic data for the Mississauga, ON, Canada common garden location
- Row 100 contains bioclimatic data for the Montpellier, France common garden location
- Row 101 contains bioclimatic data for the Uppsala, Sweden common garden location

2. In order to calculate climatic distance, we used the resultant PCA loadings for each T. repens population of origin and each common garden location. Climatic distance was calculated using only the first three dimensions, because they cumulatively account for 97% of variation, but all six dimensions are included in the dataset for completeness. Albano_et_al_Climatic_Distance_PCA_Loadings.csv

The six dimensions are labeled as Dim.1 through Dim.6 in columns A through F. Again, rows are not labeled because the subsequent code required for generating PCA figures will not run properly if the dataframe does not contain only numerical values and column headers.
- Rows 2 through 48 contain PCA loadings for all North American T. repens populations of origin
- Rows 49 through 97 contain PCA loadings for all European T. repens populations of origin
- Row 98 contains PCA loadings for the Lafayette, LA, USA common garden location
- Row 99 contains PCA loadings for the Mississauga, ON, Canada common garden location
- Row 100 contains PCA loadings for the Montpellier, France common garden location
- Row 101 contains PCA loadings for the Uppsala, Sweden common garden location

3. All identifying, genotyping, and fitness data for all individual T. repens plants in each common garden can be found in Albano_et_al_Experiment_Data.csv The
columns contained within this dataset are explained below.

- Garden: The location of the common garden (Toronto, Lafayette, Uppsala, or Montpellier) in which each T. repens individual was planted
- Position: The unique position of each individual T. repens plant within each common garden location
- Population: The identification number for each sampled T. repens population. Populations 001 through 050 are from North America, ordered from furthest south to furthest north. Populations 056 through 159 are scattered throughout Europe.
- Rep: Replicate number within each population within each common garden location (1 through 5 in Toronto and Montpellier, 1 through 8 in Lafayette and Uppsala)
- ContOrigin: The continent from which each T. repens individual originated (NAM or EUR)
- ContGarden: The continent on which the common garden is found into which each T. repens individual was planted
- Latitude: Latitude of origin of each population. Units = decimal degrees
- Longitude: Longitude of origin of each population. Units = decimal degrees
- ChangeLat: Change in latitude for each plant between its location of origin and the common garden location. Units = decimal degrees
- Cyanotype: Describes 4 T. repens phenotypes (AcLi, Acli, acLi, or acli) with respect to HCN and its underlying loci based on the presence/absence of a dominant allele at each of two loci (Ac/ac or Li/li)
- Cyanogenesis: Describes whether or not a plant is able to produce HCN (AcLi = HCN+, Acli/acLi/acli = HCN-)
- Ac_ac: Describes whether or not a plant contains at least one dominant allele at the Ac/ac locus (Ac = yes, ac = no)
- Li_li: Describes whether or not a plant contains at least one dominant allele at the Li/li locus (Li = yes, li = no)
- Survived1st: Whether or not each individual plant survived the first growing season. Survived = 1, Did not survive = 0. *Note: the 1st growing season was 2020 for NAM and 2021 for EUR*
- MaxArea1st: The maximum lateral area taken up by each plant throughout the monthly measurements of the first growing season. Units = cm^2
- GrowthRate1st: Growth rate of each plant during the first growing season, measured as ln(maximum lateral area) minus the ln(initial lateral area) divided by the # of days between those timepoints. Units = cm^2/day
- Herbivory1st: The % of leaf area consumed by herbivores, averaged across 5 trifoliate leaves per plant
- FloweredFull: Whether or not each individual plant flowered during the full experiment (both growing seasons). Flowered = 1, Did not flower = 0.
- SeededFull: Whether or not each individual plant produced seed during the full experiment (both growing seasons). Produced seed = 1, Did not produce seed = 0.
- FlowerHeadNumberFull: The total number of flower heads produced by each individual plant throughout the full experiment (both growing seasons)
- SeedSetMassFull: The total seed set mass produced by the flower heads collected from each individual plant throughout the full experiment (both growing seasons). Units = grams
- WeightedCD: The climatic distance between each population of origin and each common garden location, calculated using the 3 most influential dimensions in a principal component analysis containing 6 bioclimatic variables of interest and weighted by the % of variation explained by each of the 3 most influential dimensions.
- MATdiff: The difference in mean annual temperature (bio1, not z-standardized) between each T. repens population of origin and each common garden location
- bio1diff: The difference in mean annual temperature (WorldClim bio1, z-standardized) between each T. repens population of origin and each common garden location
- bio5diff: The difference in mean temperature of the warmest month (WorldClim bio5, z-standardized) between each T. repens population of origin and each common garden location
- bio6diff: The difference in mean temperature of the warmest month (WorldClim bio6, z-standardized) between each T. repens population of origin and each common garden location
- bio12diff: The difference in annual precipitation (WorldClim bio12, z-standardized) between each T. repens population of origin and each common garden location
- bio13diff: The difference in precipitation of the wettest month (WorldClim bio13, z-standardized) between each T. repens population of origin and each common garden location
- bio14diff: The difference in precipitation of the driest month (WorldClim bio14, z-standardized) between each T. repens population of origin and each common garden location
- MA_1st_z: (z-standardized) The maximum lateral area taken up by each plant throughout the monthly measurements of the first growing season
- GR_1st_z: (z-standardized) Growth rate of each plant during the first growing season, measured as ln(maximum lateral area) minus the ln(initial lateral area) divided by the # of days between those timepoints
- FH_Full_z: (z-standardized) The total number of flower heads produced by each individual plant throughout the full experiment (both growing seasons)
- SS_Full_z: (z-standardized) The total seed set mass produced by the flower heads collected from each individual plant throughout the full experiment (both growing seasons)
