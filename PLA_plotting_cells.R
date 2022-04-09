rm(list = ls()) # clearing the global work environment

#Load packages:
library(readr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(psych)

#Read in .csv output file from CellProfiler:
#replace /path/to/folder/ below with path to experiment folder
dir_path <- "/path/to/folder/"
expt_name <- "EMD-LMNB1_MCF7" # name of experiment folder   
root_dir <- paste(dir_path, expt_name, "/", sep="")
mydata <- read_csv(paste(root_dir, expt_name, "_Cells.csv", sep=""))   


#Add columns with information on slide and well number and experimental conditions based on lab book notes 
#Create a new column "Slide" with information on slide number from image file name
mydata$Slide <- ifelse(mydata$ImageNumber <=13 , c("S1"), c("S2"))
mydata$Slide[mydata$ImageNumber >=25] <-"S3"

#Create a new column "Well" with well names A1-3 and B1-3 from image file name.
mydata$Well <- ifelse(mydata$ImageNumber <=3 , c("A2"), c("A3"))
mydata$Well[mydata$ImageNumber >=7] <-"B2"
mydata$Well[mydata$ImageNumber >=10] <-"B3"
mydata$Well[mydata$ImageNumber >=13] <-"A2"
mydata$Well[mydata$ImageNumber >=16] <-"A3"
mydata$Well[mydata$ImageNumber >=19] <-"B2"
mydata$Well[mydata$ImageNumber >=22] <-"B3"
mydata$Well[mydata$ImageNumber >=25] <-"A2"
mydata$Well[mydata$ImageNumber >=28] <-"A3"
mydata$Well[mydata$ImageNumber >=31] <-"B2"
mydata$Well[mydata$ImageNumber >=34] <-"B3"

#Create a new column "Condition" with information on experimental condition
mydata$Condition[mydata$Slide =="S1"] <-"EMD only"
mydata$Condition[mydata$Slide == "S2"] <-"LMNB1 only"
mydata$Condition[mydata$Slide == "S3"] <-"EMD-LMNB1 complex"


#Change "Condition" to a factor instead of a character and define the levels to get it in the right order in the plots later
mydata$Condition <-factor(mydata$Condition, levels = c("EMD only", "LMNB1 only", "EMD-LMNB1 complex"))
class(mydata$Condition)
levels(mydata$Condition)

#Rename column
mydata <- mydata %>%
  rename(BlobCount = Children_BlobsInCells_Count)


#Summary of statistics
summary_stat <- describeBy(mydata$BlobCount, mydata$Condition, quant=c(0.25,0.75), mat=TRUE)
summary_stat


#Make Kruskal Wallis test
kruskal <-kruskal.test(BlobCount ~ Condition, data = mydata) 
kruskal

#Make Dunn's test
stat.test <- mydata %>%
  dunn_test(BlobCount ~ Condition, p.adjust.method = "bonferroni")%>%
  add_significance()
stat.test

#Export descriptive statistics and significance to .csv:
write.table(summary_stat, file = "EMD_LMNB1_MCF7_PLA_sum_stat.csv", sep = ",", row.names = FALSE )
write.table(stat.test, file = "EMD_LMNB1_MCF7_PLA_Signif.csv", sep = ",", row.names = FALSE )


#Rename Condition for the plot
mydata$Condition <- recode_factor(mydata$Condition, 'EMD'='Protein A only', 'LMNB1'='Protein B only', 'EMD-LMNB1'='Protein A-Protein B complex' )


#split data set for plotting based on condition
EMD_LMNB1<-mydata %>%
  filter(Condition == "EMD-LMNB1 complex")

EMD_LMNB1_omitting<-mydata %>%
  filter(Condition != "EMD-LMNB1 complex")


#Plot RCP/cell
ymax <- max(mydata$BlobCount) # ymax for plots

Plot1<-ggplot(EMD_LMNB1, aes(x=Condition,y=BlobCount)) +
  geom_boxplot(aes(fill=Condition)) + 
  ylim(0, ymax)+
  theme_classic()+
  scale_fill_grey(start=0.9, end=0.5)+
  labs(fill = "", x= "", y="RCPs/cell", title="PLA EMD + LMNB1 complex") 
Plot1

Plot2<-ggplot(EMD_LMNB1_omitting, aes(x=Condition,y=BlobCount)) +
  geom_boxplot(aes(fill=Condition)) + 
  ylim(0, ymax)+
  theme_classic()+
  scale_fill_grey(start=0.9, end=0.5)+
  labs(fill = "", x= "", y="RCPs/cell", title="PLA EMD + LMNB1 omitting controls") 
Plot2


#Export plots to pdf
ggsave("EMD_LMNB1_complex_PLA_Plot.pdf", plot = Plot1, width=4, height=5 )
ggsave("EMD_LMNB1_omitting_PLA_Plot.pdf", plot = Plot2, width=4, height=5 )

