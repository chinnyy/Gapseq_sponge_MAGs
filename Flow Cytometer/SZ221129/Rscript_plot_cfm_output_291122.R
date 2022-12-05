# Load data 

# setwd
setwd("C:/Users/chiny/Desktop/UNSW/CELL_COUNT/Flow Cytometer/SZ221129")

library(tidyr)
library(dplyr) # data wrangling
library(ggplot2) # visualization
library(scales) # visualization

avg_data = read.csv(file= "Flow_cytometer_output_for_R_AVG_291122.csv", header = TRUE, sep = ",", dec = ".", fill = TRUE, fileEncoding = 'UTF-8-BOM')

avg_data<- avg_data[,1:7]

avg_long <- gather(avg_data, stain, cell_count, STY09:DS_Live, factor_key=TRUE)

# Plot for the full set of data

# Regular version
ggplot (avg_long,
        aes(x=CFU,
            y=cell_count))+
  geom_line(aes(colour=stain))+
  geom_point(aes(colour=stain))+
  theme_light()

# log version 
ggplot (avg_long,
        aes(x=log10(CFU),
            y=log10(cell_count)))+
  geom_line(aes(colour=stain))+
  geom_point(aes(colour=stain))+
  ylab(bquote(Cell~Count/mL (Log[10])))+
  xlab(bquote(CFU/mL (Log[10])))+
  #ylim(0,10)+
  #xlim(0,10)+
  scale_y_continuous(breaks=seq(0,12,2))+
  scale_x_continuous(breaks=seq(0,12,2))+
  ylim(0,10)+
  xlim(0,10)+
  theme_light()

# Cut of -inf

avg_data_ss <- subset(avg_long, log_CFU != "-Inf")

# log version 
ggplot (avg_data_ss,
        aes(x=log10(CFU),
            y=log10(cell_count)))+
  geom_line(aes(colour=stain))+
  geom_point(aes(colour=stain))+
  ylab(bquote(Cell~Count/mL (Log[10])))+
  xlab(bquote(CFU/mL (Log[10])))+
  #ylim(0,10)+
  #xlim(0,10)+
  scale_y_continuous(breaks=seq(0,12,2))+
  scale_x_continuous(breaks=seq(0,12,2))+
  ylim(0,10)+
  xlim(0,10)+
  theme_light()

# Logged data

avg_long$log_CFU <- log10(avg_long[,c("CFU")])  

avg_long$log_cc <- log10(avg_long[,c("cell_count")]) 

# Plot for a smaller dataset removing CFU 10^6-10^8

# subsetting the dataset

avg_data_ss <- avg_data[4:10,]
                      
avg_long_ss <- gather(avg_data_ss, stain, cell_count, STY09:CFU_SZ, factor_key=TRUE)

avg_long_ss$log_CFU <- log10(avg_long_ss[,c("CFU")])  

avg_long_ss$log_cc <- log10(avg_long_ss[,c("cell_count")]) 

# Plot for the subset of data
ggplot (avg_long_ss,
        aes(x=CFU,
            y=cell_count))+
  geom_line(aes(colour=stain))+
  geom_point(aes(colour=stain))+
  theme_light()

ggplot (avg_long_ss,
        aes(x=log10(CFU),
            y=log10(cell_count)))+
  geom_line(aes(colour=stain))+
  geom_point(aes(colour=stain))+
  ylab(bquote(Cell~Count/mL (Log[10])))+
  xlab(bquote(CFU/mL (Log[10])))+
  theme_light()

### LINE PLOTS WITH GROUPING TO SEE ACCURACY OF DATA

raw_data = read.csv(file= "Flow_cytometer_output_for_R_291122.csv", header = TRUE, sep = ",", dec = ".", fill = TRUE, fileEncoding = 'UTF-8-BOM')

raw_long <- gather(raw_data, stain, cell_count, STY09:DS, factor_key=TRUE)

ggplot(raw_long, aes(x=as.factor(CFU), y=log10(cell_count))) + 
  geom_boxplot(aes(colour=stain))+
  facet_wrap(~stain)+
  ylab(bquote(Cell~Count/mL (Log[10])))+
  xlab(bquote(CFU/mL))+
  theme_light()+
  theme(legend.position="none")
