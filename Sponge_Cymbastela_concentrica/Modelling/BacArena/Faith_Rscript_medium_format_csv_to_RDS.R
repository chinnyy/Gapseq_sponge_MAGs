##### medium formatting Rscript ##########
### Purpose is to change a regular csv output from medium prediction into a formatted version and convert it to a .RDS file

# Set working directory
setwd <- setwd("C:/Users/chiny/Desktop/UNSW/Gapseq_sponge_MAGs/Sponge_Cymbastela_concentrica/Modelling/qsub_gapseq") # Change this accordingly to what sponge you are working on

bac_id<- 'bin.2' # Change this accordingly to what the file where the medium csv is stored in

model_id<- 'bin.2_mod_2-medium' # Change this accordingly to what the file where the medium csv is stored in

# Load the medium .csv file
medium<- read.csv(paste(setwd,'/',bac_id,'/',model_id,'.csv',sep = ''))

# Renaming the columns 
colnames(medium) <- c("Exchange","Description","Input_mM")

# Create dataframe containing additional information

df <- data.frame (Input_mol  = rep(c(""),times=nrow(medium)),
                  init_conc = rep(c("\"NA\""),times=nrow(medium)),
                  init_area  = rep(c(5),times=nrow(medium)),
                  D_fec = rep(c(0),times=nrow(medium)),
                  D_muc  = rep(c(""),times=nrow(medium)),
                  D_omu = rep(c(""),times=nrow(medium)),
                  D_imu  = rep(c(""),times=nrow(medium)),
                  D_inter = rep(c(""),times=nrow(medium)),
                  D_hum  = rep(c(""),times=nrow(medium)),
                  A_fec = rep(c(""),times=nrow(medium)),
                  A_omu  = rep(c(""),times=nrow(medium)),
                  A_imu = rep(c(""),times=nrow(medium)),
                  A_hum = rep(c(""),times=nrow(medium)),
                  A_inter  = rep(c(""),times=nrow(medium)),
                  pde = rep(c(""),times=nrow(medium))
                  )
# Append it to the medium to create the formatted version 
medium_formatted<- cbind(medium,df)

# Saving the new formatted medium as a CSV
write.csv(medium_formatted,paste(bac_id,'/',model_id,'-medium_formatted.csv',sep = ''), row.names = F)

# Making a BacArena medium from Gapseq medium and saving it:
saveRDS(medium_formatted, file = paste(bac_id,'/',model_id,'-medium_formatted.RDS',sep = ''))

