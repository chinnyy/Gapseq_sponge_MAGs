# Get formatted reaction equations and fluxes  by providing MS ID (e.g.rxnxxxx_c0): 
#### IMPORTANT: this script can be applied to any input file with ModelSEED id, e.g. rxnxxxx_c0.
#### IMPORTANT: ModelSEED reaction id must be in the 1st column with header start with 'react_id'.

# Set to where you want the script to run on 
dir_out = setwd('/Users/chiny/Desktop/UNSW/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Modelling/Scripts')

gapseq_version  = 'gapseq_v20220329'

# Import the gap-filled modified object model (gfM): 
# Change accordingly to the name of the object model!
gfM_dir = 'Actino_5_Burgsdorf_2021_full/'
model_id = 'Actino_5_Burgsdorf_2021_full' 

# Navigate to the python code "Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py"
python_code  = 'py /Users/chiny/Desktop/UNSW/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Modelling/Scripts/Python_scripts/Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py'

# Navigate to where the seed reaction tsv "seed_reactions_corrected.formatted_reat_add.tsv"
dict_seed_reactions_corrected_tsv   = paste('/Users/chiny/Bioinfo/software/gapseq/',gapseq_version,'/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv',sep = '')

# Navigate to where the Smat_fluxes file is located (e.g."050_model_XXXX_gfM_Smat_fluxes.csv")
infile_flux = paste(dir_out,'/07_gf_model/',gfM_dir,'050_model_',model_id,'_gfM_Smat_fluxes.csv',sep = '')

# Navigate to file with flux of MS ids (050_model_xxxxxx_gfMOM_Smat_fluxes.csv)
infile_key_list = paste(dir_out,'/07_gf_model/',gfM_dir,model_id,'_completely_filled_in_reactions.txt',sep = '')

outfile_key_list_equation_formatted = paste(dir_out,'/medium',model_id,'_completely_filled_in_reactions.equation_formatted_n_fluxes.txt',sep = '')

commands = paste(python_code,
                 '-infile_MSid',infile_key_list,
                 '-infile_flux', infile_flux, 
                 '-db',dict_seed_reactions_corrected_tsv,
                 '-outfile',outfile_key_list_equation_formatted,
                 '-gs',gapseq_version, sep = ' ')

system(commands)


