rm(list = ls()) # clearing the global work environment

library("readr")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("rstatix")
library("psych")

# reading in csv file produced by 'MB_thresholding_cells.ipynb'
# wherein all blobs have been categorised as one of the two antibodies or both
# replace /path/to/folder/ below with path to experiment folder
dir_path <- "/path/to/folder/"
expt_name <- "EMD-LMNB1_MCF7" # name of experiment folder
root_dir <- paste(dir_path, expt_name, "/", sep="")
mydata <- read_csv(paste(root_dir, expt_name, "_IntegratedIntensity.csv", sep=""))

# adding Categorization as a factor made from Category
mydata$Categorization[mydata$Category =="Both"] <- "Interaction"
mydata$Categorization[mydata$Category == "TexRed"] <- "Free EMD"
mydata$Categorization[mydata$Category == "Atto647"] <- "Free LMNB1"
# putting levels in the right order for plotting
mydata$Categorization <-factor(mydata$Categorization, levels = c("Interaction", "Free EMD", "Free LMNB1"))

title_plot1 <- "EMD-LMNB1"
title_plot2 <- "EMD only control"
title_plot3 <- "LMNB1 only control"
ymax <- max(mydata$BlobCount) # ymax for plots
legend_position <- c(0.8, 0.9)

class(mydata$Categorization)
levels(mydata$Categorization)

mydataMouse <- mydata %>%
  filter(Condition == "EMD")

mydataRabbit <- mydata %>%
  filter(Condition == "LMNB1")

mydataBoth <- mydata %>%
  filter(Condition == "EMD-LMNB1")

# make plot1 - boxplots (both antibodies)
plot<-ggplot(mydataBoth, aes(x=Categorization,y=BlobCount)) +
  geom_boxplot(aes(fill=Categorization), show.legend = TRUE) + 
  theme_classic()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position=legend_position, plot.title = element_text(hjust = 0.5))+
  scale_fill_grey(start=0.3, end=0.9)+
  ylim(0, ymax)+
  labs(fill = "", x= "", y="RCPs/cell", title=title_plot1) 

plot
save_name <- paste(root_dir, title_plot1, ".pdf", sep="")
ggsave(save_name, device="pdf", width=4, height=5)

# make plot2 - boxplots (Mouse antibody)
plot<-ggplot(mydataMouse, aes(x=Categorization,y=BlobCount)) +
  geom_boxplot(aes(fill=Categorization), show.legend = TRUE) + 
  theme_classic()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position=legend_position, plot.title = element_text(hjust = 0.5))+
  scale_fill_grey(start=0.3, end=0.9)+
  ylim(0, ymax)+
  labs(fill = "", x= "", y="RCPs/cell", title=title_plot2) 

plot
save_name <- paste(root_dir, title_plot2, ".pdf", sep="")
ggsave(save_name, device="pdf", width=4, height=5)

# make plot3 - boxplots (Rabbit antibody)
plot<-ggplot(mydataRabbit, aes(x=Categorization,y=BlobCount)) +
  geom_boxplot(aes(fill=Categorization), show.legend = TRUE) + 
  theme_classic()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position=legend_position, plot.title = element_text(hjust = 0.5))+
  scale_fill_grey(start=0.3, end=0.9)+
  ylim(0, ymax)+
  labs(fill = "", x= "", y="RCPs/cell", title=title_plot3) 

plot
save_name <- paste(root_dir, title_plot3, ".pdf", sep="")
ggsave(save_name, device="pdf", width=4, height=5)

# compute and save summary statistics
x <- describeBy(mydataBoth$BlobCount, mydataBoth$Categorization, quant=c(0.25,0.75), mat=TRUE)
x
save_name <- paste(root_dir, title_plot1, ".csv", sep="")
write.table(x, file=save_name, sep=",", row.names=FALSE)

x <- describeBy(mydataMouse$BlobCount, mydataMouse$Categorization, quant=c(0.25,0.75), mat=TRUE)
x
save_name <- paste(root_dir, title_plot2, ".csv", sep="")
write.table(x, file=save_name, sep=",", row.names=FALSE)

x <- describeBy(mydataRabbit$BlobCount, mydataRabbit$Categorization, quant=c(0.25,0.75), mat=TRUE)
x
save_name <- paste(root_dir, title_plot3, ".csv", sep="")
write.table(x, file=save_name, sep=",", row.names=FALSE)

