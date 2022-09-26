################################################################ katana ####################################################################################################
############ 
library(lattice)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)
library(sybil)
library(Matrix)
library(BacArena)
library(parallel)
library(plyr)

# ############################ Download RDS file and run on MacOS ############################
setwd <- setwd()

# Change blow:
dir_group = '/Actino_mtf_220602_full_recipe_iter_15/'

# import BacArena RDS file:
simulation_loop <- readRDS(paste(getwd,dir_group,'/BacArena_STY_Actino_400grids_mineral_sw_rl_Actino.RDS',sep = ''))
# ############################################################################################
# 08. get mflux of each species 
# How to choose the time point?
# check the 001_plotGrowthCurve_biomass_STY_Actino_mineral_sw_rl.png, 
# and find the time point when the curve has the largest slope.
timept = 7


# Change blow:
# Extract all the rxn IDs(= ModelSEED ID = MS ID) of the model.
write.table(ldply(simulation_loop[[1]]@mfluxlist[[timept]][["gapseq version: 1.2 5c13d4b; Sequence DB md5sum: 1139b8e (2022-04-01, Bacteria)"]],data.frame), file = paste(getwd,dir_group,'/080_simulation_mfluxlist_tp',timept,'_STY_Merged_Actino.tsv',sep = ''), quote = F,sep = '\t',row.names = F)

# 08.1New. Get formatted reaction equations and fluxes in BacArena results by providing MS ID (e.g.rxnxxxx_c0): 
#### IMPORTANT: this script can be applied to any input file with ModelSEED id, e.g. rxnxxxx_c0.
#### IMPORTANT: ModelSEED reaction id must be in the 1st column with header start with 'react_id', while the 2nd column include fluxes.
#### Updates on 05/08/21: provide a file "seed_reactions_corrected.formatted_reat_add.status.??.tsv".(e.g. seed_reactions_corrected.formatted_reat_add.status.Actino.tsv (a file that was copied from the directory of the gap-filled modified model (gfMOM), and add the OTU name to this file))
#### Updates on 07/08/21: for the "Single species" analysis, provide directory of dictionary and directory of output files separately.  
# e.g. python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_formatted_equations_n_flux_of_MS_react_in_gfModel_v2.py -in_dir -dict_dir -timept -dict_cus_react -gs_version -model_id
# python_code         = '/Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_formatted_equations_n_flux_of_MS_react_in_gfModel_v3_BacArena.py'
python_code         = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/07_Python_scripts/Get_formatted_equations_n_flux_of_MS_react_in_gfModel_v3_BacArena.py'

#Change here:

workdate = '20210825'

input_dir       = paste(getwd, dir_group, sep = '')
current_dir     = getwd()
dict_dir        = paste(current_dir,'/',sep = '')
dict_cus_react  = paste('/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/06_DB_files/Reaction_id_in_ReactionPool_Addreact_',workdate,'.txt',sep = '')
gs_version      = 'v20210504' #gapseq_version
model_id        = 'Actino'

commands = paste('python3',python_code,'-in_dir',input_dir, '-dict_dir', dict_dir,'-tp', timept, '-dict',dict_cus_react,'-gs',gs_version,'-modNm',model_id, sep = ' ')
system(commands)


# 08.2. Get formatted reaction equations and fluxes by providing MS ID (e.g.rxnxxxx_c0): 
# e.g. python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_EX_cpd_annotate_BacArena.py -in_dir -timept -gs_version -model_id
python_code         = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/07_Python_scripts/Get_EX_cpd_annotate_BacArena.py'
#Change here:
input_dir       = paste(getwd, dir_group, sep = '')
gs_version      = 'v20220329' #gapseq_version
model_id        = 'Actino'

commands = paste('python3',python_code,'-in_dir',input_dir,'-tp', timept, '-gs',gs_version,'-modNm',model_id, sep = ' ')
system(commands)


######################## Follow with python script of 'Get_ETC_complexes_reaction_flux_argv.py '######################## 
# Get_active_rxn_Equ_annot_by_met_from_final_model_argv.py
### update @2021-07-01: 
#### Give results of flux values of ETC chain complexes reactions
##### Input key file: /Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/EC_num_of_ETC_complexes.txt 
###### could be updates if reactions were missing, or EC num is wrong.
#### update @2021-08-09:
###### Input file of Bacarena result.

# E.g.
# python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_ETC_complexes_reaction_flux_BacArena_argv.py -in_dict -in_key -gs
mycmd1 = 'python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_ETC_complexes_reaction_flux_BacArena_argv.py'
mycmd2 = paste('-in_dict ', getwd, dir_group, '/080_simulation_mfluxlist_tp', timept, '_STY_Merged_Actino.Equ2_lt0.tsv', sep = '')
mycmd5 = paste('-in_key ', '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/06_DB_files/EC_num_of_ETC_complexes.txt',sep = '')

system(paste(mycmd1,mycmd2,mycmd5,sep = ' '))
