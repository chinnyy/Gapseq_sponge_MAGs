# @2022-05-24
# run /Users/zzfanyi/PycharmProjects/EMP500/BacArena/Get_active_gapfilled_rxn_in_BacArena_argv.py

######## 1. For the Actino_mtf_220602_full_recipe_iter_15:
cmds = 'python3 /Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/07_Python_scripts/Get_active_gapfilled_rxn_in_BacArena_argv.py'

flag1 = '-mydir_infile_BA'
cmd1 = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/05_BacArena/Single_species_test/Actino_mtf_220602_full_recipe_iter_15/'

flag2 = '-infile_BA'
cmd2 = '080_simulation_mfluxlist_tp7_STY_Merged_Actino.Equ2_lt0.tsv'

flag3 = '-mydir_infile_gs'
cmd3 = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/03_Get_ObjectModel/'

flag4 = '-infile_gs'
cmd4 = '20220602_medium_Actino_5_Burgsdorf_2021_pos_completely_filled_in_reactions.equation_formatted_n_fluxes.txt'

flag5 = '-outfile'
cmd5 = '080_simulation_mfluxlist_tp7_STY_Merged_Actino.Equ2_lt0.gapfilled.tsv'

system(paste(cmds,flag1,cmd1,flag2,cmd2,flag3,cmd3,flag4,cmd4,flag5,cmd5,sep = ' '))

# You will get an output file "080_simulation_mfluxlist_tp7_STY_Merged_Actino.Equ2_lt0.gapfilled.tsv",
# which includes all reactions that were filled to the gapseq model and activated in the single speceis 
# simulation using BacArena.


