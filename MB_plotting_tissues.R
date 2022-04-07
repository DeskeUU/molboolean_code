rm(list = ls()) # clearing the global work environment

library("readr")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("rstatix")
library("psych")

# reading in csv file produced by 'MB_thresholding_tissues.ipynb'
# wherein all blobs have been categorised as one of the two antibodies or both
# replace /path/to/folder/ below with path to experiment folder
dir_path <- "/path/to/folder/"
expt_name <- "Ecad-bcat_prostate" # name of experiment folder
root_dir <- paste(dir_path, expt_name, "/", sep="")
mydata <- read_csv(paste(root_dir, expt_name, "_IntegratedIntensity.csv", sep=""))

# adding Categorization as a factor made from Category
mydata$Categorization[mydata$Category =="Both"] <- "Interaction"
mydata$Categorization[mydata$Category == "TexRed"] <- "Free E-cad"
mydata$Categorization[mydata$Category == "Atto647"] <- "Free b-cat"

# putting levels in the right order for plotting
mydata$Categorization <-factor(mydata$Categorization, levels = c("Interaction", "Free E-cad", "Free b-cat"))

title_plot1 <- "E-cad-b-cat"
title_plot2 <- "E-cad only control"
title_plot3 <- "b-cat only control"

class(mydata$Categorization)
levels(mydata$Categorization)

mydataMouse <- mydata %>%
  filter(Condition == "Ecad")

mydataRabbit <- mydata %>%
  filter(Condition == "bcat")

mydataBoth <- mydata %>%
  filter(Condition == "Ecad-bcat")

# make plot1 - piechart (both antibodies)
df = mydataBoth[, c("BlobCount", "Categorization")]
df2 = aggregate(.~Categorization,data=df,FUN=sum)

df3 <- df2
df3$Percentage <- df3$BlobCount / sum(df3$BlobCount) * 100
save_name <- paste(root_dir, title_plot1, ".csv", sep="")
write.table(df3, file=save_name, sep=",", row.names=FALSE)

plot<-ggplot(df2, aes(x='',y=BlobCount, fill=Categorization)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_grey(start=0.3, end=0.9)+
  theme(axis.line=element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))+
  labs(title=title_plot1)

plot
save_name <- paste(root_dir, title_plot1, ".pdf", sep="")
ggsave(save_name, device="pdf", width=4, height=3)

# make plot2 - piechart (Mouse antibody)
df = mydataMouse[, c("BlobCount", "Categorization")]
df2 = aggregate(.~Categorization,data=df,FUN=sum)

df3 <- df2
df3$Percentage <- df3$BlobCount / sum(df3$BlobCount) * 100
save_name <- paste(root_dir, title_plot2, ".csv", sep="")
write.table(df3, file=save_name, sep=",", row.names=FALSE)

plot<-ggplot(df2, aes(x='',y=BlobCount, fill=Categorization)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_grey(start=0.3, end=0.9)+
  theme(axis.line=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(title=title_plot2)

plot
save_name <- paste(root_dir, title_plot2, ".pdf", sep="")
ggsave(save_name, device="pdf", width=4, height=3)

# make plot3 - piechart (Rabbit antibody)
df = mydataRabbit[, c("BlobCount", "Categorization")]
df2 = aggregate(.~Categorization,data=df,FUN=sum)

df3 <- df2
df3$Percentage <- df3$BlobCount / sum(df3$BlobCount) * 100
save_name <- paste(root_dir, title_plot3, ".csv", sep="")
write.table(df3, file=save_name, sep=",", row.names=FALSE)

plot<-ggplot(df2, aes(x='',y=BlobCount, fill=Categorization)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_grey(start=0.3, end=0.9)+
  theme(axis.line=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(title=title_plot3)

plot
save_name <- paste(root_dir, title_plot3, ".pdf", sep="")
ggsave(save_name, device="pdf", width=4, height=3)

