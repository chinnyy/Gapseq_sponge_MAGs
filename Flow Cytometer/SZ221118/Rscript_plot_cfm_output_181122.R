# Load data 

# setwd
setwd("C:/Users/chiny/Desktop/UNSW/CELL_COUNT/Flow Cytometer/SZ221118")

library(tidyr)
library(dplyr) # data wrangling
library(ggplot2) # visualization

avg_data = read.csv(file= "Flow_cytometer_output_for_R__AVG_181122.csv", header = TRUE, sep = ",", dec = ".", fill = TRUE, fileEncoding = 'UTF-8-BOM')

avg_long <- gather(avg_data, stain, cell_count, STY09:CFU_SZ, factor_key=TRUE)

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
  theme_light()

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


