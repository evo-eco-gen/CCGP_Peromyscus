#### Plot age of repeat elemnts based on ParseRM.pl 'age mode' results
## Prepare the workspace:
rm(list=ls())
gc()
setwd("/Users/kozakk/Documents/CCGP/Peromyscus/Genome/Repeatable/ages")

## read in the ages table estimated with multiple rounds of ParseRM.pl
## perl parseRM.pl -i pman.all_masked.align -a $age -m 0.0022 -v
ages <- read.csv("P.m.sonoriensis.ParseRM.repeatAges.csv", header=T)
head(ages)
colnames(ages)[5] <- "percentage_masked"
plot(ages$Age, ages$percentage_masked)

## plot properly with ggplot2
library(tidyverse)
ages_data <- ggplot(ages, aes(x = Age, y = percentage_masked))
ages_data + geom_point(shape = 23, fill="dark grey") + theme_light() + xlab("Myr ago") + ylab("% genome (no redundancy)")

## Save
ggsave("repeat_ages.pdf", plot = last_plot(), dpi=300)
save.image(file="plot_repeat_ages.RData")

## Tidy up
cat("\014")
q()
