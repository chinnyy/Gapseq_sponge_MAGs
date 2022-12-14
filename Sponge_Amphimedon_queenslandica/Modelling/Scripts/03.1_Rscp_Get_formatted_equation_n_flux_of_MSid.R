# Get formatted reaction equations and fluxes  by providing MS ID (e.g.rxnxxxx_c0): 
#### IMPORTANT: this script can be applied to any input file with ModelSEED id, e.g. rxnxxxx_c0.
#### IMPORTANT: ModelSEED reaction id must be in the 1st column with header start with 'react_id'.

# python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py -infile_MSid -outfile -infile_flux -gs -db
dir_out = getwd()

# gapseq_version      = 'gapseq_v20210504'
gapseq_version      = 'gapseq_v20220329'

python_code         = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/Python_scripts/Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py'
dict_seed_reactions_corrected_tsv   = paste('/Users/zzfanyi/Bioinfo/software/gapseq/',gapseq_version,'/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv',sep = '')

# Change this:
infile_flux = '/Users/zzfanyi/Bioinfo/PostDoc/Gordon_and_Betty_Moore_Foundation/Petrosia_ficiformis/Laura_DTE_techniques/2022-04-16_Actinobacteria/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/07_gf_model/20220602_medium_Actino_5_Burgsdorf_2021_pos/050_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes.csv'

infile_key_list                     = paste(dir_out,'/20220602_medium_Actino_5_Burgsdorf_2021_pos_completely_filled_in_reactions.txt',sep = '')
outfile_key_list_equation_formatted = paste(dir_out,'/20220602_medium_Actino_5_Burgsdorf_2021_pos_completely_filled_in_reactions.equation_formatted_n_fluxes.txt',sep = '')

commands = paste('python3',python_code,'-infile_MSid',infile_key_list,'-infile_flux', infile_flux, '-db',dict_seed_reactions_corrected_tsv,'-outfile',outfile_key_list_equation_formatted,'-gs',gapseq_version, sep = ' ')
system(commands)
