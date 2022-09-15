# @2021-07-06.
#### IMPORTANT: this script can be applied to any input file with ModelSEED id, e.g. rxnxxxx_c0.
#### IMPORTANT: ModelSEED reaction id must be in the 1st column with header start with 'react_id'.


import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_formatted_equations_n_flux_of_MS_react_in_gfModel/Get_formatted_equations_n_flux_of_MS_react_in_gfModel.py -infile_MSid -outfile -infile_flux -gs -db')
parser.add_argument('-infile_MSid', required=True, help='file with MS ids provided in the 1st column (rxnxxxx_c0)')
parser.add_argument('-outfile', required=True, help='output file, ends with .equation_formatted_n_fluxes.txt')
parser.add_argument('-infile_flux', required=True, help='file with flux of MS ids (050_model_xxxxxx_gfMOM_Smat_fluxes.csv)')

parser.add_argument('-gs', required=False, default='gapseq_v20210504', help='gapseq version, default. gapseq_v20210504')
parser.add_argument('-db', required=False, default='/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_v20210504/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv', help='path to seed_reactions_corrected.tsv')
args = vars(parser.parse_args())

infile_key_list = args['infile_MSid']
outfile_key_list_equation = args['outfile']
infile_flux = args['infile_flux']
gapseq_version = args['gs']
dict_seed_reactions_corrected_tsv = args['db']

##################### Or run script within pycharm:
######################### user define start #########################
#
# gapseq_version                      = 'gapseq_v20210504'
# # input 1: dict file of 'seed_reactions_corrected.formatted_reat_add.tsv'. e.g '/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_v20210504/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'
# dict_seed_reactions_corrected_tsv   = '/Users/zzfanyi/Bioinfo/software/gapseq/'+gapseq_version+'/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'
# # input 2: file with MS ids provided in the 1st column (rxnxxxx_c0):
# dir_infile                          = '/Users/zzfanyi/Bioinfo/Sponge_Petrosia/To_Shan_2021-03-03_from_Ilia/selected_by_Shan/Archaea_2/Step_protocals_using_R/03_Get_ObjectModel/20210704_diet_sw_Archaea_2_arc/'
# infile_key_list                     = dir_infile + '20210704_diet_sw_Archaea_2_arc_completely_filled_in_reactions.txt'
# outfile_key_list_equation           = dir_infile + '20210704_diet_sw_Archaea_2_arc_completely_filled_in_reactions.equation_formatted_n_fluxes.txt'
# # input 3: file with flux of MS ids (050_model_xxxxxx_gfMOM_Smat_fluxes.csv).
# dir_infile_flux = '/Users/zzfanyi/Bioinfo/Sponge_Petrosia/To_Shan_2021-03-03_from_Ilia/selected_by_Shan/Archaea_2/Step_protocals_using_R/07_gf_model/20210704_diet_sw_Archaea_2_arc/'
# infile_flux = dir_infile_flux + '050_model_Archaea2_gfMOM_Smat_fluxes.csv'


# gapseq_version                      = 'gapseq_v20210504'
# # input 1: dict file of 'seed_reactions_corrected.formatted_reat_add.tsv'. e.g '/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_v20210504/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'
# dict_seed_reactions_corrected_tsv   = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Workflow_in_silico_medium_prediction/7_Dict_files/seed_reactions_corrected.formatted_reat_add.tsv'
# # input 2: file with MS ids provided in the 1st column (rxnxxxx_c0):
# dir_infile                          = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Workflow_in_silico_medium_prediction/2_Output_files/My_Input_Archaea_2_react_diff/'
# infile_key_list                     = dir_infile + '20210704_diet_sw_Archaea_2_arc_completely_filled_in_reactions.txt'
# outfile_key_list_equation           = dir_infile + '20210704_diet_sw_Archaea_2_arc_completely_filled_in_reactions.equation_formatted_n_fluxes.txt'
# # input 3: file with flux of MS ids (050_model_xxxxxx_gfMOM_Smat_fluxes.csv).
# dir_infile_flux = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Workflow_in_silico_medium_prediction/2_Output_files/My_Input_Archaea_2-draft/'
# infile_flux = dir_infile_flux + '050_model_Archaea_2-draft_model_Smat_fluxes.csv'
#
#
# gapseq_version                      = 'gapseq_v20210504'
# # input 1: dict file of 'seed_reactions_corrected.formatted_reat_add.tsv'. e.g '/Users/zzfanyi/Bioinfo/software/gapseq/gapseq_v20210504/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'
# dict_seed_reactions_corrected_tsv   = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/gapseq_v20220329/Shan_output_files/seed_reactions_corrected.formatted_reat_add.tsv'
# # input 2: file with MS ids provided in the 1st column (rxnxxxx_c0):
# dir_infile                          = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/03_Get_ObjectModel/'
# infile_key_list                     = dir_infile + '20220602_medium_Actino_5_Burgsdorf_2021_pos_completely_filled_in_reactions.txt'
# outfile_key_list_equation           = dir_infile + '20220602_medium_Actino_5_Burgsdorf_2021_pos_completely_filled_in_reactions.equation_formatted_n_fluxes.txt'
# # input 3: file with flux of MS ids (050_model_xxxxxx_gfMOM_Smat_fluxes.csv).
# dir_infile_flux = '/Users/zzfanyi/Gh_Sflabelliformis_8_MAGs/Gapseq_sponge_MAGs/Sponge_Petrosia_ficiformis/Actino_5_Burgsdorf_2021/03_Gapseq_models/Actino_5_Burgsdorf_2021/Gram_pos/2022-04-04_Model_checkup_pkg/07_gf_model/20220602_medium_Actino_5_Burgsdorf_2021_pos/'
# infile_flux = dir_infile_flux + '050_model_Actino_5_Burgsdorf_2021_gfM_Smat_fluxes.csv'


######################## user define end #########################

# st0. create a seed reaction flux dict using '050_model_xxxxxx_gfMOM_Smat_fluxes.csv'
dict_seed_rxn_id_2_fluxes = {}
for each in open(infile_flux):
    each_split = each.strip().split(',')
    # print(each_split)
    if each.startswith(','):
        header_dict_0 = each_split
        # print(header_dict_0)

    else:
        # print(each)
        rxn_id_dict = each_split[0].split('_')[0]
        rxn_flux = each_split[1]

        # rxn_name_dict = each.strip().split('\t')[2]
        # rxn_equation_dict = each.strip().split('\t')[6]
        # rxn_equation_def_dict = each.strip().split('\t')[7]
        # rxn_description_dict = [rxn_id_dict,rxn_name_dict,rxn_equation_dict,rxn_equation_def_dict]
        # print(rxn_description_dict)
        # rxn_description_dict = dict_rxn_id_2_description['\t'.join(rxn_id_dict)]
        dict_seed_rxn_id_2_fluxes[rxn_id_dict] = rxn_flux

# print(len(dict_seed_rxn_id_2_description))
print('%s%s%s%s' % ('st0. There are ',len(dict_seed_rxn_id_2_fluxes),' MS IDs with associated fluxes found in the file ', infile_flux))
# st0. There are 1668 MS IDs with associated fluxes found in the file /Users/zzfanyi/Bioinfo/Sponge_Petrosia/To_Shan_2021-03-03_from_Ilia/selected_by_Shan/Archaea_2/Step_protocals_using_R/07_gf_model/20210704_diet_sw_Archaea_2_arc/050_model_Archaea2_gfMOM_Smat_fluxes.csv



# st1. create a seed reaction dictionary using dat/seed_reactions_corrected.tsv
dict_seed_rxn_id_2_description = {}
for each in open(dict_seed_reactions_corrected_tsv):
    each_split = each.strip().split('\t')
    if each.startswith('Rid\t'):
        header_dict = each_split
    else:
        # print(each)
        rxn_id_dict = each_split[0].split('_Shan')[0]
        rxn_description_dict = each_split

        # rxn_name_dict = each.strip().split('\t')[2]
        # rxn_equation_dict = each.strip().split('\t')[6]
        # rxn_equation_def_dict = each.strip().split('\t')[7]
        # rxn_description_dict = [rxn_id_dict,rxn_name_dict,rxn_equation_dict,rxn_equation_def_dict]
        # print(rxn_description_dict)
        # rxn_description_dict = dict_rxn_id_2_description['\t'.join(rxn_id_dict)]
        dict_seed_rxn_id_2_description[rxn_id_dict] = '\t'.join(rxn_description_dict)

# print(len(dict_seed_rxn_id_2_description))
print('%s%s%s%s' % ('st1. There are ',len(dict_seed_rxn_id_2_description),' MS IDs with associated equations found in gapseq ', gapseq_version))
# st1. There are 34822 MS IDs with associated equations found in gapseq v20210318

# print(dict_seed_rxn_id_2_description)



#st2. Find if interested MS ids (from list) are in OM.
output_handle = open(outfile_key_list_equation, 'w')
for each in open(infile_key_list):
    each_split = each.strip().split('\t')
    # print(each_split)
    if each.startswith('react_id'):
        header = each_split
        for_print = '%s\t%s\t%s\n' % ('\t'.join(header),'\t'.join(header_dict),'sol$fluxes')
        # print(for_print)
        output_handle.write(for_print)
    else:
        MS_id = each_split[0].split('_')[0]
        # print(MS_id)
        react_equation = 'MS id not found in the gapseq database'
        react_flux = 'MS id not found in the final model'
        if MS_id in dict_seed_rxn_id_2_description:
            react_equation = dict_seed_rxn_id_2_description[MS_id]
            rxn_flux = dict_seed_rxn_id_2_fluxes[MS_id]

        for_print = ('%s\t%s\t%s\n' % ('\t'.join(each_split), react_equation, rxn_flux))
        output_handle.write(for_print)
output_handle.close()

