## SET UP
# Load relevant packages
library(dplyr)

# Set working directory
wd = setwd('/Users/chiny/Desktop/UNSW/Gapseq_sponge_MAGs/Sponge_Cymbastela_concentrica/Modelling/Scripts/07_gf_model')


##############

# Step 1: Determine what compound is in model 1 but not in model 2 and vice versa

# In this case I will be comparing between bin.1_full and bin.1_mod_1, please change the name of the models accordingly to where it is stored or names

bin_dir = "/bin.1/"

mod_1_name = "bin.1_full"
  
mod_2_name = "bin.1_mod_1"


# Load the two files "binXXX-medium.csv" that you want to compare between

mod_1_medium <- read.csv(paste(wd,bin_dir,mod_1_name,"/",mod_1_name,"-medium.csv",sep = ""), header=TRUE)

mod_2_medium <- read.csv(paste(wd,bin_dir,mod_2_name,"/",mod_2_name,"-medium.csv",sep = ""), header=TRUE)

# Compare the differences between compounds in 2 media
# This combines the compounds that are found in model 1 but not in model 2 and vice versa

diff_med <- rbind(anti_join(mod_1_medium,mod_2_medium, by= "compounds"),anti_join(mod_2_medium,mod_1_medium, by= "compounds"))

##############

# Extract relevant reactions according to the different compounds 
# This is done by looking into the column "met.id" and seeing if the compound is there

# Load the 2 relevant "xxxx_completely_filled_in_reactions.equation_formatted_n_fluxes.txt" involved

mod_1_rxn <- read.delim(paste(wd,bin_dir,mod_1_name,"/",mod_1_name,"_completely_filled_in_reactions.equation_formatted_n_fluxes.txt",sep = ""))
  
mod_2_rxn <- read.delim(paste(wd,bin_dir,mod_2_name,"/",mod_2_name,"_completely_filled_in_reactions.equation_formatted_n_fluxes.txt",sep = ""))

# Change the compound name in diff_med to match that in mod_1_rxn

diff_med$c0_cpd<- paste(diff_med[,1],'[c0]',sep = "")


# Create new folder where the output will be stored
out_dir<- paste(wd,bin_dir,"Compare_",mod_1_name,'_and_',mod_2_name,sep = "")
dir.create(out_dir)

for (i in 1:nrow(diff_med)){
  # For mod 1
  mod_1_rxn$index <- grepl(diff_med[i,4], mod_1_rxn$met_id, fixed=TRUE) # Get the index of reactions with the compounds mentioned
  filtered_1 <- mod_1_rxn[mod_1_rxn$index == "TRUE",] # Keep the rows with "TRUE"
  write.csv(filtered_1,paste(out_dir,"/",mod_1_name,'_completely_filled_in_reactions.equation_formatted_n_fluxes.',diff_med[i,1],'_',diff_med[i,2],'.csv',sep = ""), row.names=FALSE)
  
  # For mod 2
  mod_2_rxn$index <- grepl(diff_med[i,4], mod_2_rxn$met_id, fixed=TRUE) # Get the index of reactions with the compounds mentioned
  filtered_2 <- mod_2_rxn[mod_2_rxn$index == "TRUE",] # Keep the rows with "TRUE"
  write.csv(filtered_2,paste(out_dir,"/",mod_2_name,'_completely_filled_in_reactions.equation_formatted_n_fluxes.',diff_med[i,1],'_',diff_med[i,2],'.csv',sep = ""), row.names=FALSE)
}








