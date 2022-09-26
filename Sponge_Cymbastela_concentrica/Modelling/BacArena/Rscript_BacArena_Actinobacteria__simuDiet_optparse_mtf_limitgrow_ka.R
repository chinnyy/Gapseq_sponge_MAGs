################################################################# katana ####################################################################################################
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
######################## argument ########################
option_list = list(
  optparse::make_option(c("-o", "--otu_name"),                 type="character",      default='Actino',       help="OTU_name, defalt: Actino"),
  optparse::make_option(c("-m", "--infile_my_model"),          type="character",      default='/srv/scratch/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/20210515_STY_Merged_OTU06_mineral_sw_3_NoH2S/STY_Merged_OTU06.RDS',            
                        help="directory of input final model, defalt: /srv/scratch/z5095298/sponge_modeling/GapSeq_v20210429_STY8taxa_gapseq_v0429_MKSdb_Diet_sw_lALL_b_uncon/OTU06_renamed/20210515_STY_Merged_OTU06_mineral_sw_3_NoH2S/STY_Merged_OTU06.RDS"),
  optparse::make_option(c("-D", "--infile_BacarenaRDS_diet"),  type="character",      default="/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3only.RDS",
                        help="Diet RDS file for bacarena with full path, defalt: /srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified_9999_NH3only.RDS"),
  optparse::make_option(c("-k", "--keywd"),                    type="character",      default="sulfite",    help="Nutirent keyword used in the folder name of output files, defalt: sulfite"),
  
  # optparse::make_option(c("-n", "--no2_con"),   type="double",       default=0.5,       help="NH3 concentration, defalt: 0.5"),
  optparse::make_option(c("-a", "--arena_mn"),  type="double",       default=20,        help="integer indicating the length of an arena, defalt: 20"),
  optparse::make_option(c("-d", "--death_r"),   type="double",       default=0,         help="A percentage of biomass reduce due to the nutrient limitation, defalt: 0"),
  optparse::make_option(c("-i", "--inocc_no"),  type="double",       default=0.05,      help="inocculum, defalt: 0.05"),
  optparse::make_option(c("-c", "--cl_no"),     type="double",       default=4,         help="Number of replicates, defalt: 4"),
  optparse::make_option(c("-t", "--iter"),      type="double",       default=200,       help="Number of iteration, defalt:200"));

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

otu_name        = opt$otu_name
infile_my_model = opt$infile_my_model
infile_diet     = opt$infile_BacarenaRDS_diet
keywd           = opt$keywd
grid_no   = opt$arena_mn
death_r   = opt$death_r
Inocc_no  = opt$inocc_no
Cl_no     = opt$cl_no
iter_no   = opt$iter

#####################################################################################################################################################################

############ Katana ############ 
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic.RDS')
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_full.RDS')
# diet <- readRDS('/srv/scratch/z5095298/software/gapseq/Diets_Sponges_July_2019/Diet_Sponges_July_2019_dl20190807_Shan_v20210425/Diet_Sponges_July_2019_dl20190807_Shan_v20210425_brief_v2_IDcheck_mineral_sw_210514/Diet_for_BacArena/mineral_sw_210514_OTU08_3_BacArena_realistic_simplified.RDS')
diet <- readRDS(infile_diet)
#####################################################################################################################################################################
getwd <- getwd()
getwd

# create a new folder:
new_folder = paste(otu_name,'_',keywd,'_iter_',iter_no,'/',sep = '')
dir.create(paste(getwd,'/',new_folder,sep = ''))


# NOB:
my_model <- readRDS(infile_my_model) #d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__UBA8639;g__bin75;s__

replicates <- Cl_no
cores <- Cl_no
cl <- makeCluster(cores, type="PSOCK")
clusterExport(cl, c("diet","my_model","getwd","new_folder","grid_no","death_r","Inocc_no","iter_no","replicates","cores"))
clusterEvalQ(cl, sink(paste0(getwd,'/',new_folder, Sys.getpid(), ".txt")))

simlist <- parLapply(cl, 1:replicates, function(i){
  sybil::SYBIL_SETTINGS("SOLVER", "cplexAPI") #https://github.com/euba/BacArena/issues/152
  print("====================================================================================")
  print("============================= Actinobacteria =======================================")
  print("====================================================================================")
  bacterium_Actino <- BacArena::Bac(my_model,lyse = F,deathrate=death_r,maxweight = 1, minweight=0.125,cellweight_mean = 0.5,cellweight_sd = 0.25, growtype='exponential',setAllExInf=FALSE,limit_growth = T)
  arena <- BacArena::Arena(n=grid_no,m=grid_no)
  inocc <- arena@n * arena@m * Inocc_no # 5% inocculum
  arena <- BacArena::addOrg(object = arena, specI = bacterium_Actino, amount = inocc*1) # amount can be added by fraction.
  arena <- BacArena::addSubs(object = arena, mediac = diet$Exchange, smax = diet$Input_mM, unit = "mM", add = T)
  arena@tstep <- 1
  simulation <- BacArena::simEnv(object = arena, time = iter_no, sec_obj = "mtf", continue = T) #minimize total flux
  # simulation <- BacArena::simEnv(object = arena, time = iter_no, continue = T) #FBA; iter 1


})
stopCluster(cl)

saveRDS(simlist, file = paste(getwd,'/',new_folder,'/BacArena_STY_',otu_name,'_400grids_mineral_sw_rl_Actino.RDS',sep = ''))

# ############################ Download RDS file and run on MacOS ############################
getwd <- getwd()
getwd
simulation_loop <- readRDS(paste(getwd,'/',new_folder,'/BacArena_STY_',otu_name,'_400grids_mineral_sw_rl_Actino.RDS',sep = ''))
# simulation_loop <- readRDS(paste(getwd,'/Actino_mtf_220602_full_recipe_iter_15/BacArena_STY_Actino_400grids_mineral_sw_rl_Actino.RDS',sep = ''))
# plotGrowthCurve(simulation_loop,use_biomass = T)
# simulation_loop <- readRDS(paste(getwd,'/Actino_mtf_220602_Thiamin_sulfo_iter_15/BacArena_STY_Actino_400grids_mineral_sw_rl_Actino.RDS',sep = ''))
# plotGrowthCurve(simulation_loop,use_biomass = T)

# # 00. plotGrowthCurve:
# # Plot growth curve for several simulation_loops.The function plotGrowthCurve takes a list of simulation_loops and plots the time course of species with standard deviation.
pdf(paste(getwd,'/',new_folder,'/001_plotGrowthCurve_biomass_STY_',otu_name,'_mineral_sw_rl.pdf',sep = ''),width = 5, height = 4)
plotGrowthCurve(simulation_loop,use_biomass = T)
dev.off()
use_biomass <- plotGrowthCurve(simulation_loop,use_biomass = T,ret_data = T)
write.csv(use_biomass, file = paste(getwd,'/',new_folder,'/001_plotGrowthCurve_biomass_STY_',otu_name,'_mineral_sw_rl.csv',sep = ''))
#
pdf(paste(getwd,'/',new_folder,'/002_plotGrowthCurve_number_STY_',otu_name,'_mineral_sw_rl.pdf',sep = ''),width = 5, height = 4)
plotGrowthCurve(simulation_loop,use_biomass = F)
dev.off()
Nouse_biomass <- plotGrowthCurve(simulation_loop,use_biomass = F,ret_data = T)
write.csv(Nouse_biomass, file = paste(getwd,'/',new_folder,'/002_plotGrowthCurve_number_STY_',otu_name,'_mineral_sw_rl.csv',sep = ''))
# To get better plot, run 081: BacArena_STY_5taxa_400grids_spongeDiet_6iterations.R


######################################## exclude biomass ########################################
####time0
type_xy_md_time0 <- evalArena(simulation_loop, plot_items = "Population", phencol = F, retdata = T, 
                              time = (seq_along(simulation_loop@simlist) - 1),show_legend = F)
write.table(type_xy_md_time0$Population$time0, file = paste(getwd,'/',new_folder,'/081_ggplot_',auto_num,'.tsv',sep = ''))

dt_type_xy_md_time0 <- as.data.frame(read.table(file = paste(getwd,'/',new_folder,'/081_ggplot_',auto_num,'.tsv',sep = '')))
# head(dt_type_xy_md_time0)
# str(dt_type_xy_md_time0)

library('ggplot2')
pdf(paste(getwd,'/',new_folder,'/081_ggplot_',auto_num,'.pdf',sep = ''),width = 5, height = 4)
ggplot(dt_type_xy_md_time0, aes(x=dt_type_xy_md_time0$x, y=dt_type_xy_md_time0$y)) + 
  geom_point(color=dt_type_xy_md_time0$type,size=0.3)+
  labs(title = "STY_8species_time0") + 
  coord_fixed(ratio=1) + 
  theme(plot.background = element_blank())+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position=c(1,1), legend.justification=c(1,1))
dev.off()
