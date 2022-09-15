# st6.4. Extract all reactions from the gap-filled model (gfM).

##################### st1. Get input model files and output docs ##################### 
library(stringr)
library(sybil)

getwd = getwd()
getwd

# Import the gap-filled modified object model (gfM): 
gfM_dir = '20220602_medium_Actino_5_Burgsdorf_2021_pos/' # manually change this!

Actino_5_Burgsdorf_2021_gfM = paste(getwd,'/Input_files/',gfM_dir,'Actino_5_Burgsdorf_2021.RDS',sep = '')

# create a new folder:
new_folder_nm = paste('/07_gf_model/',gfM_dir,sep = '')
dir.create(paste(getwd,new_folder_nm,sep = ''))

# Assign output RDS file name:
outfile_add_043= paste(getwd, new_folder_nm,'043_model_Actino_5_Burgsdorf_2021_gfM_list_react.tsv',sep = '')
outfile_add_045= paste(getwd, new_folder_nm,'045_model_Actino_5_Burgsdorf_2021_gfM_list_metabolites.tsv',sep = '')
outfile_add_046= paste(getwd, new_folder_nm,'046_model_Actino_5_Burgsdorf_2021_gfM_Smat.csv',sep = '')
infile_046     = outfile_add_046
outfile_050    = paste(getwd, new_folder_nm,'050_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes.csv',sep = "")
outfile_050_2  = paste(getwd, new_folder_nm,'050_2_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes.csv',sep = "")
outfile_050_3 = paste(getwd, new_folder_nm,'050_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes_ETC_complexes.tsv',sep = "")

gapseq_version = 'gapseq_v20220329'
model_id = 'Actino_5_Burgsdorf_2021'

model_Actino_5_Burgsdorf_2021_gfM = readRDS(Actino_5_Burgsdorf_2021_gfM)
model_Actino_5_Burgsdorf_2021_gfM
# model name:             Actino_5_Burgsdorf_2021 
# number of compartments  3 
# c0 
# e0 
# p0 
# number of reactions:    2829 
# number of metabolites:  2409 
# number of unique genes: 1505 
# objective function:     +1 EX_cpd11416_c0
########################## st4. Extract all reactions from the modified model (MM). #############################  
# 043.What are those reactions in the model?
model_Actino_5_Burgsdorf_2021_gfM@react_num
# [1] 1976

# Extract all reactions from the addReact OM.
model_Actino_5_Burgsdorf_2021_gfM_list <- list(model_Actino_5_Burgsdorf_2021_gfM@react_id, model_Actino_5_Burgsdorf_2021_gfM@react_name, model_Actino_5_Burgsdorf_2021_gfM@react_rev, model_Actino_5_Burgsdorf_2021_gfM@react_single, model_Actino_5_Burgsdorf_2021_gfM@react_de,
                                               model_Actino_5_Burgsdorf_2021_gfM@lowbnd, model_Actino_5_Burgsdorf_2021_gfM@uppbnd, model_Actino_5_Burgsdorf_2021_gfM@obj_coef, model_Actino_5_Burgsdorf_2021_gfM@gprRules, model_Actino_5_Burgsdorf_2021_gfM@gpr,
                                               model_Actino_5_Burgsdorf_2021_gfM@react_attr)
df_model_Actino_5_Burgsdorf_2021_gfM_list  <- as.data.frame(model_Actino_5_Burgsdorf_2021_gfM_list, row.names = NULL)
colnames(df_model_Actino_5_Burgsdorf_2021_gfM_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr'
                                                         ,'seed',"rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart"
                                                         ,"send","pathway","blast.status","pathway.status","complex","exception","complex.status","gs.origin","annotation","MNX_ID","seedID","keggID","biggID","biocycID"
)
write.table(df_model_Actino_5_Burgsdorf_2021_gfM_list, file = outfile_add_043,quote = F,row.names = F,sep = '\t')
# The output contains many KEGG reaction ids, which can be converted to SEED RXN ids by using seed_reactions.tsv.


########################## st5. Extract all metabolites from the modified model (MM). #############################  
# 045.What are the metabolites in this new model?

# number of reactions of ******model_Actino_5_Burgsdorf_2021_gfM****: 
model_Actino_5_Burgsdorf_2021_gfM@met_num
# [1] 1899
model_Actino_5_Burgsdorf_2021_gfM_met_name = model_Actino_5_Burgsdorf_2021_gfM@met_name
model_Actino_5_Burgsdorf_2021_gfM_met_id = model_Actino_5_Burgsdorf_2021_gfM@met_id
# model_Actino_5_Burgsdorf_2021_gfM: Generate a dataframe with metabolite info for each model
model_Actino_5_Burgsdorf_2021_gfM_list_met <- list(model_Actino_5_Burgsdorf_2021_gfM@met_id, model_Actino_5_Burgsdorf_2021_gfM@met_name, model_Actino_5_Burgsdorf_2021_gfM@met_comp, model_Actino_5_Burgsdorf_2021_gfM@met_single, model_Actino_5_Burgsdorf_2021_gfM@met_de)
df_model_Actino_5_Burgsdorf_2021_gfM_list_met <- as.data.frame(model_Actino_5_Burgsdorf_2021_gfM_list_met, row.names = NULL)
colnames(df_model_Actino_5_Burgsdorf_2021_gfM_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_Actino_5_Burgsdorf_2021_gfM_list_met_sort <- df_model_Actino_5_Burgsdorf_2021_gfM_list_met[order(df_model_Actino_5_Burgsdorf_2021_gfM_list_met$met_id),]
write.table(df_model_Actino_5_Burgsdorf_2021_gfM_list_met_sort, file = outfile_add_045, quote = F,row.names = F,sep = '\t')


########################## st6. Extract Stoichiometric matrix from the modified model (MM). #############################  
# 046.What are the stiochiometric matrix (model_xxx@S) in this model?
nrow(model_Actino_5_Burgsdorf_2021_gfM@S)
# [1] 1899

ncol(model_Actino_5_Burgsdorf_2021_gfM@S)
# [1] 1976


model_Actino_5_Burgsdorf_2021_gfM_Smat <- model_Actino_5_Burgsdorf_2021_gfM@S
model_Actino_5_Burgsdorf_2021_gfM_Smat_mx = as.matrix(model_Actino_5_Burgsdorf_2021_gfM@S)
class(model_Actino_5_Burgsdorf_2021_gfM_Smat_mx)
# [1] "matrix"

colnames(model_Actino_5_Burgsdorf_2021_gfM_Smat_mx)
# NULL
model_Actino_5_Burgsdorf_2021_gfM@react_num
# [1] 1976
model_Actino_5_Burgsdorf_2021_gfM@react_id
colnames(model_Actino_5_Burgsdorf_2021_gfM_Smat_mx)<-cbind(model_Actino_5_Burgsdorf_2021_gfM@react_id)

row.names(model_Actino_5_Burgsdorf_2021_gfM_Smat_mx)
# NULL
model_Actino_5_Burgsdorf_2021_gfM@met_num
# 1899
model_Actino_5_Burgsdorf_2021_gfM@met_id
row.names(model_Actino_5_Burgsdorf_2021_gfM_Smat_mx)<-cbind(model_Actino_5_Burgsdorf_2021_gfM@met_id)

write.csv(model_Actino_5_Burgsdorf_2021_gfM_Smat_mx, file = outfile_add_046, quote = FALSE,row.names = T)
# In this table, you can check the flux (values) of a compound (e.g. cpd16503[e0]) of this model based on the reaction ID (e.g. EX_cpd16503_e0).


########################## st6. Check flux from the gfM. #############################  
############ 050. What are the active reactions which use or produce a certain compound in ******model_Actino_5_Burgsdorf_2021_gfM****:
mod = model_Actino_5_Burgsdorf_2021_gfM
dim(mod)
# [1] 1899 1976
# this will give a new column called 'sol$fluxes' with the same length to the 'reaction id', which shows fluxes of each reaction.
?optimizeProb
# more functions please read the manual of R package sybil.
dim(mod)
# [1] 2712 3147

# sol <- sybil::optimizeProb(mod, retOptSol=F, algorithm = "mtf", mtfobj = T) # find a list (of 7,sol), which include flux info of 765 reactions for this model.
sol <- sybil::optimizeProb(mod, retOptSol=F, algorithm = "fba") # find a list (of 7,sol), which include flux info of 765 reactions for this model.
model_Actino_5_Burgsdorf_2021_gfM_flux_list <- as.list(sol$fluxes)
length(model_Actino_5_Burgsdorf_2021_gfM_flux_list)
# [1] 3147

model_Actino_5_Burgsdorf_2021_gfM_flux_list <- as.data.frame(sol$fluxes)
dim(model_Actino_5_Burgsdorf_2021_gfM_flux_list)
# [1] 3147    1


model_Actino_5_Burgsdorf_2021_gfM_Smat <- mod@S

model_Actino_5_Burgsdorf_2021_gfM_Smat <- read.csv(file = infile_046,row.names = 'X')
t_model_Actino_5_Burgsdorf_2021_gfM_Smat <- as.data.frame(t(model_Actino_5_Burgsdorf_2021_gfM_Smat))
dim(t_model_Actino_5_Burgsdorf_2021_gfM_Smat)
# [1] 3147 2712

cbind_t_model_Actino_5_Burgsdorf_2021_gfM_Smat<- as.data.frame(cbind(model_Actino_5_Burgsdorf_2021_gfM_flux_list, t_model_Actino_5_Burgsdorf_2021_gfM_Smat))
dim(cbind_t_model_Actino_5_Burgsdorf_2021_gfM_Smat)
# [1] 1976 1900

write.csv(cbind_t_model_Actino_5_Burgsdorf_2021_gfM_Smat, quote = FALSE,file = outfile_050)
# IMPORTANT: check the flux value of ADD_BIOMASS, if it is negative, it is not right!

######################## Follow with python script of 'Assign_names_2_ID_of_active_rxn_in_Smax_arg.py'######################## 
############ 051. assign reaction name and metabolite name to 050 file in ******model_Actino_5_Burgsdorf_2021_gfM****:
# New on 2021-mar-30
# Assign_names_2_ID_of_active_rxn_in_Smax_arg.py
#### IMPORTANT-2: this script used data in gapseq_v20200407, should be updated with the pipline version used for modeling.
#### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)

#### update @2022-04-11: Modify the python script under /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Check_gfM/Assign_names_2_ID_of_active_rxn_in_Smax_arg.py
# IMPORTANT!!!!! Get the db file Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv by running:
# /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Shan_output_files/Get_db_files_in_Shan_output_files.py
# Run this in python (usually this only run once for a gapseq version):
# 'st4.Format_seed_reactions_Willis.py'
######################  run script in command line.
getwd
input_dir = paste(getwd,'/07_gf_model/',gfM_dir,sep = '')
model_id = 'Actino_5_Burgsdorf_2021'
gapseq_version = gapseq_version

mycmd1 = 'python3 /Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/Python_scripts/Assign_names_2_ID_of_active_rxn_in_Smax_arg.py'
# Please open "Assign_names_2_ID_of_active_rxn_in_Smax_arg.py" and modify the full path of Dict input files from line 78 to line 80.                  
system(paste(mycmd1,
             '-in_dir',input_dir,
             '-modNm', model_id, 
             '-gs',gapseq_version, sep = ' '))
# output files:
# 051_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes.flux_lt0.csv
# 052_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes.flux_lt0.annotate.tsv

################################################################################################
# Run the following script once to get the "seed_reactions_corrected.formatted_reat_add.status.tsv" for the related 043 output file:
# 'Get_a_file_of_formatted_MS_reactios_with_assigned_reaction_and_pathway_status.py'
mycmd1 = 'python3 /Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/Python_scripts/Get_a_file_of_formatted_MS_reactions_with_assigned_reaction_and_pathway_status.py'
# Please open the python script "Get_a_file_of_formatted_MS_reactions_with_assigned_reaction_and_pathway_status.py", and change the path of "dict_gs_seed_rxn" in line 63.
mycmd2 = paste('-dir_gfM ',getwd,'/07_gf_model/',gfM_dir, sep = '')
mycmd3 = '-gs gapseq_v20220329'
mycmd4 = '-dict_pwy_tbl /Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Actino_5_Burgsdorf_2021-all-Pathways.tbl'
mycmd5 = '-dict_043 043_model_Actino_5_Burgsdorf_2021_gfM_list_react.tsv'


system(paste(mycmd1,mycmd2,mycmd3,mycmd4,mycmd5,sep = ' '))
######################## Follow with python script of 'Get_active_rxn_Equ_annot_by_met_from_final_model_argv.py'######################## 
# Get_active_rxn_Equ_annot_by_met_from_final_model_argv.py
### update @2021-02-24: reaction direction in the output file is the same in file dat/seed_reactions_corrected.tsv
#### The true direction should combine with flux values (depending on either positive or negative)
##### The flux here does not equal to the fluxes in the gap-filling (depends on algorithms) step!!!
###### Py script updated on 2021-apr-14.
####### run 'st4.Format_seed_reactions_Willis.py' before running this:

# Create a list of compounds that you have interests in their flux.
list_met = list('cpd00013', #NH3-c0	and NH3-e0
                # 'cpd00165', #Hydroxylamine
                # 'cpd00011', #CO2
                # 'cpd00025', #H2O2
                # 'cpd00242', #H2CO3
                'cpd00418', # NO
                'cpd00528', # N2
                # 'cpd00073', # Urea
                'cpd00075', # Nitrite
                'cpd00209', # Nitrate
                # 'cpd00136', # 4-Hydroxybenzoate
                
                # 'cpd00020', # Pyruvate for glycolysis
                # 'cpd00027', # D-Glucose
                # 'cpd00159', # L-Lactate
                # 'cpd00221', # D-Lactate
                # 'cpd00086', # Propionyl-CoA for pyruvate fermentation to propanoate I
                # 'cpd19008', # (R)-Acetoin for pyruvate fermentation to (R)-acetoin I
                # 'cpd11640', # H2 hydrogen
                
                # 'cpd00793', # B1 Thiamine phosphate
                # 'cpd00220', #	B2 (Riboflavin); Riboflavin-c0
                # 'cpd00644', # B5 (pantothenate); PAN, Pantothenat
                # 'cpd00263', #	B6 (Pyridoxine)
                # 'cpd00104', #	B7 (Biotin); BIOT
                # 'cpd00730', #	B12 OH-Cbl (Hydroxocobalamin)
                # 'cpd01826', #	B12 CN-Cbl (Cyanocobalamin)
                # 'cpd12878', #	B12 Me-Cbl (Methylcobalamin)
                # 'cpd00166', #	B12 Ado-Cbl (adenosylcobalamin); Calomide-c0
                # 'cpd00218', # Niacin
                # 'cpd00133', # Nicotinamide
                #
                # 'cpd03177', #"Hydrazine"
                # 'cpd11798', #Ferrocytochrome c2
                # 'cpd01024', #methane
                # 'cpd00146', # Carbamoylphosphate[0]: product of the 1st reaction of the urea cycle
                # 'cpd10516', # fe3+
                'cpd00007', #O2
                # 'cpd00002', #ATP
                # 'cpd00067', # H+
                # 'cpd11606', #Menaquinone 7
                # 'cpd15560' #Ubiquinone-8
                'cpd00074', # Sulfer      S (0)
                'cpd00239', # H2S         H2S (-1)
                'cpd00081', # Sulfite     H2SO3 (+4)
                'cpd00268', # Thiosulfate H2S2O3 (+2)
                'cpd00048', # Sulfate     H2SO4 (+3)
                # 'cpd01270', # FMNH2
                'cpd00193', # APS:Adenosine 5'-phosphosulfate
                'cpd00042', # GSH: 5-L-Glutamyl-L-cysteinylglycine
                'cpd17456', # S-Sulfanylglutathione, 4
                'cpd01414', # Tetrathionate
                'cpd00109', # CytC 3+
                'cpd00210',  # Taurine      organic S
                'cpd11416' #biomass
                
)


mycmd1 = 'python3 /Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/Python_scripts/Get_active_rxn_Equ_annot_by_met_from_final_model_argv3.py'
mycmd2 = paste('-in_dir ',getwd,'/07_gf_model/',gfM_dir,sep = '')
mycmd6 = paste('-gs ',gapseq_version,sep = '')
mycmd7 = paste('-modNm ', model_id,sep='')

for (met in list_met) {
  print(met)
  mycmd4 = paste('-met ',met,sep = '')
  system(paste(mycmd1,mycmd2,mycmd4,mycmd6,mycmd7,sep = ' '))
  
}
######################## Follow with python script of 'Get_ETC_complexes_reaction_flux_argv.py '######################## 
# Get_active_rxn_Equ_annot_by_met_from_final_model_argv.py
### update @2021-07-01: 
#### Give results of flux values of ETC chain complexes reactions
##### Input key file: /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt 
###### could be updates if reactions were missing, or EC num is wrong.
# E.g.
# python3 /Users/zzfanyi/PycharmProjects/my_pys/EMP500/GapSeq/Get_ETC_complexes_reaction_flux_argv.py 
# -in_dict /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfM_Smat_fluxes.csv 
# -in_key /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt 
# -gs gapseq_v20210409 
# -out /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/STY_models/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU07_renamed/Step_protocals_using_R/07_gf_model/20210630_diet_sw_no_oxygen_STY_Merged_OTU07-draft__12_thiosulfate/050_model_OTU07_gfM_Smat_fluxes_ETC_complexes.tsv'

mycmd1 = 'python3 /Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/Python_scripts/Get_ETC_complexes_reaction_flux_argv.py'
mycmd2 = paste('-in_dict ',outfile_050, sep = '')
# mycmd3 = '-file 052_model_OTU07_gfM_Smat_fluxes.flux_lt0.annotate.tsv'
mycmd5 = paste('-in_key ','/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/06_DB_files/EC_num_of_ETC_complexes.txt',sep = '')
mycmd6 = paste('-gs ',gapseq_version,sep = '')
mycmd7 = paste('-out ', outfile_050_3, sep='')

system(paste(mycmd1,mycmd2,mycmd5,mycmd6,mycmd7,sep = ' '))