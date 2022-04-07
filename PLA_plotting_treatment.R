rm(list = ls()) # clearing the global work environment


#Load packages:
library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(psych)


#Read in .csv output file from CellProfiler:
# wherein all blobs have been categorised as one of the two antibodies or both
# replace /path/to/folder/ below with path to experiment folder
dir_path <- "/path/to/folder/"
expt_name <- "Ecad-bcat_TGFb" # name of experiment folder  
root_dir <- paste(dir_path, expt_name, "/", sep="")
mydata <- read_csv(paste(root_dir, expt_name, "_Cells.csv", sep="")) 

#Add columns with information on slide and well number and experimental conditions based on lab book notes 
#Create a new column "Slide" with information on slide number from image file name
mydata$Slide[mydata$ImageNumber == 1 ] <-"S1"
mydata$Slide[mydata$ImageNumber == 2 ] <-"S2"
mydata$Slide[mydata$ImageNumber == 3 ] <-"S3"
mydata$Slide[mydata$ImageNumber == 4 ] <-"S1"
mydata$Slide[mydata$ImageNumber == 5 ] <-"S2"
mydata$Slide[mydata$ImageNumber == 6 ] <-"S3"
mydata$Slide[mydata$ImageNumber == 7 ] <-"S1"
mydata$Slide[mydata$ImageNumber == 8 ] <-"S2"
mydata$Slide[mydata$ImageNumber == 9 ] <-"S3"
mydata$Slide[mydata$ImageNumber == 10 ] <-"S1"
mydata$Slide[mydata$ImageNumber == 11 ] <-"S2"
mydata$Slide[mydata$ImageNumber == 12 ] <-"S3"
mydata$Slide[mydata$ImageNumber == 13 ] <-"S1"
mydata$Slide[mydata$ImageNumber == 14 ] <-"S2"
mydata$Slide[mydata$ImageNumber == 15 ] <-"S3"
mydata$Slide[mydata$ImageNumber == 16 ] <-"S1"
mydata$Slide[mydata$ImageNumber == 17 ] <-"S2"
mydata$Slide[mydata$ImageNumber == 18 ] <-"S3"

#Create a new column "Well" with well names A1-3 and B1-3 from image file name.
mydata$Well <- ifelse(mydata$ImageNumber <=3 , c("A1"), c("A2"))
mydata$Well[mydata$ImageNumber >=7] <-"A3"
mydata$Well[mydata$ImageNumber >=10] <-"B1"
mydata$Well[mydata$ImageNumber >=13] <-"B2"
mydata$Well[mydata$ImageNumber >=16] <-"B3"

#Create a new column "Condition" with information on experimental condition
mydata$Condition[mydata$ImageNumber <10] <-"Treatment"
mydata$Condition[mydata$ImageNumber >=10] <-"Control"

#Change "Well", "Slide" and "Condition" to factors instead of characters and defining the levels to get it in the right order in the plots later
mydata$Well <-factor(mydata$Well, levels = c("A1", "B1","A2", "B2", "A3", "B3"))
class(mydata$Well)
levels(mydata$Well)

mydata$Slide <-factor(mydata$Slide, levels = c("S1", "S2", "S3"))
class(mydata$Slide)
levels(mydata$Slide)

mydata$Condition <-factor(mydata$Condition, levels = c("Control", "Treatment"))
class(mydata$Condition)
levels(mydata$Condition)


#Rename column
mydata <- mydata %>%
  rename(BlobCount = Children_BlobsInCells_Count)


#Summary of statistics
summary_stat <- describeBy(mydata$BlobCount, mydata$Condition, quant=c(0.25,0.75), mat=TRUE)
summary_stat

#Make Wilcoxon rank sum test
stat.test <- mydata %>%
  wilcox_test(BlobCount ~ Condition) %>%
  add_significance()
stat.test 


#Export descriptive statistics and significance to .csv:
write.table(summary_stat, file = "Ecad_bcat_TGFb_Hacat_PLA_sum_stat.csv", sep = ",", row.names = FALSE )
write.table(stat.test, file = "Ecad_bcat_TGFb_Hacat_PLA_Signif.csv", sep = ",", row.names = FALSE )


#Plot RCP/cell
Plot<-ggplot(mydata, aes(x=Condition,y=BlobCount)) +
  geom_boxplot(aes(fill=Condition)) + 
  theme_classic()+
  scale_fill_grey(start=0.9, end=0.5)+
  labs(fill = "", x= "", y="RCP/cell", title="PLA Ecad + bcat TGFb" ) 

#Add p-values to plot
stat.test <- stat.test %>% add_xy_position(x = "Condition")
Plot <-Plot + stat_pvalue_manual(stat.test)


#Export plots to pdf
ggsave("Ecad_bcat_TGFb_PLA_Plot.pdf", plot = Plot, width=4, height=5 )






