rm(list = ls()) # clearing the global work environment

library("readr")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("rstatix")
library("psych")

# reading in csv file produced by 'MB_thresholding_treatments.ipynb'
# wherein all blobs have been categorised as one of the two antibodies or both
# replace /path/to/folder/ below with path to experiment folder
dir_path <- "/path/to/folder/"
expt_name <- "Ecad-bcat_TGFb" # name of experiment folder
root_dir <- paste(dir_path, expt_name, "/", sep="")
mydata <- read_csv(paste(root_dir, expt_name, "_IntegratedIntensity.csv", sep=""))

# adding Treatment as a factor made from Condition
mydata$Treatment[mydata$Condition =="Treatment"] <- "Treatment"
mydata$Treatment[mydata$Condition == "Control"] <- "Control"

# putting levels in the right order for plotting
mydata$Treatment <-factor(mydata$Treatment, levels = c("Treatment", "Control"))
title_plot <- expt_name

class(mydata$Treatment)
levels(mydata$Treatment)

mydata$Categorization[mydata$Category =="Both"] <- "Interaction"
mydata$Categorization[mydata$Category == "TexRed"] <- "Free E-cadherin"
mydata$Categorization[mydata$Category == "Atto647"] <- "Free b-catenin"

# putting levels in the right order for plotting
mydata$Categorization <-factor(mydata$Categorization, levels = c("Interaction", "Free E-cadherin", "Free b-catenin"))

class(mydata$Categorization)
levels(mydata$Categorization)

# make boxplots with significance added
p <- ggplot(mydata, aes(x=Treatment,y=BlobCount,color=Categorization, fill=Categorization)) +
  geom_boxplot(colour = 'black') +
  facet_grid(. ~ Categorization) +
  labs(fill = "", x= "", y="RCP proportion/cell", title=title_plot )

p

stat.test <- mydata %>%
  group_by(Categorization) %>%
  dunn_test(BlobCount ~ Treatment, p.adjust.method = "bonferroni")%>%
  add_significance()

stat.test
save_name <- paste(root_dir, title_plot, "_Significance.csv", sep="")
write.table(stat.test, file=save_name, sep=",", row.names=FALSE)


stat.test <- stat.test %>% add_xy_position(x = "Treatment")

p + stat_pvalue_manual(stat.test, step.group.by = "Categorization") +
  theme_classic() +
  scale_fill_grey(start=0.3, end=0.9) +
  theme(strip.text = element_blank(), plot.title = element_text(hjust = 0.5))

save_name <- paste(root_dir, title_plot, ".pdf", sep="")
ggsave(save_name, device="pdf", width=7, height=5)

# compute and save summary stats
mydataTreated <- mydata %>%
  filter(Condition == "Treatment")

x <- describeBy(mydataTreated$BlobCount, mydataTreated$Categorization, quant=c(0.25,0.75), mat=TRUE)
x

save_name <- paste(root_dir, title_plot, "_Treated.csv", sep="")
write.table(x, file=save_name, sep=",", row.names=FALSE)

mydataControl <- mydata %>%
  filter(Condition == "Control")

x <- describeBy(mydataControl$BlobCount, mydataControl$Categorization, quant=c(0.25,0.75), mat=TRUE)
x

save_name <- paste(root_dir, title_plot, "_Control.csv", sep="")
write.table(x, file=save_name, sep=",", row.names=FALSE)


