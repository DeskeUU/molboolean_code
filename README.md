# molboolean_code

In this repository we provide code to accompany the paper:

## MolBoolean: A Boolean Analysis of Protein-Protein Interactions at a Molecular Level

### By Raykova et al.

The code is based on both python and R and includes Jupyter notebooks (https://jupyter.org/) for setting the thresholds, as discussed below. The code has been tested on macOS version 12.1, python version 3.7.6 and R version 3.6.2.

Python libraries required: pandas, numpy, matplotlib.

R libraries required: readr, ggplot2, dplyr, ggpubr, rstatix, psych.

## Instructions, demonstration and reproduction

The code is split by 'cells', 'tissues' and 'treatment'. For all three cases the RCPs in the cells or tissues (termed 'blobs' in the code) are thresholded as belonging to either the Mouse, Rabbit or both antibodies. This thresholding is done based on inspecting the plots of the intensities of the stains, with Mouse (TexRed) on the x-axis and Rabbit (Atto647) on the y-axis. Also used to inform this thresholding is a histogram of the angles for the blobs in the x-y plots, with an angle of zero for a blob parallel with the x-axis and an angle of 90 for a blob parallel with the y-axis. See the Jupyter notebooks (.ipynb files) for examples for each of the three cases. For some cases a new orgin was defined for thresholding of the intensities (as shown here in MB_thresholding_treatments.ipynb).

Once the blobs have been so categorised a csv file is produced. This csv file is run through the R scripts to produce the final plots as shown in the manuscript: boxplots for the 'cells', piecharts for the 'tissues' and boxplots with signficance added for the 'treatments'.

The examples given include EMD-LMNB1 (Fig. 3a and Supp. Fig 4a in manuscript) for the 'cells', Ecad-bcat in prostate (Fig. 6c and Supp. Fig. 7c) for the 'tissues' and TGFb1 treated Ecad-bcat (Fig. 5a) for the 'treatments'.

Reproducing the results as shown in the Jupyter notebooks and the manuscript plots requires creating csv files from the three CellProfiler output files (with the suffices 'Image', 'Both_Cells' and 'BlobsInCells') from the source data shared with the manuscript.

Running the code for the thresholding and plotting in each case takes only a couple of minutes.
