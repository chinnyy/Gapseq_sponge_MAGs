
#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=30gb
#PBS -l walltime=23:59:00
#PBS -M z5434493@ad.unsw.edu.au
#PBS -m ae

#################### gapseq #################### 

# Version gapseq_v20210501 new system, new R lib directory.

module load perl/5.28.0
module load git/2.22.0
module load bedtools/2.27.1
module load blast+/2.9.0
module load hmmer/3.2.1
module load glpk/4.65
module load barrnap/0.9
module load gcc/7.3.0
module load exonerate/2.2.0
module load parallel/20190522
module load libsbml/5.18.0 # module load since 2020-12-29
# module add R/3.5.3 
module add R/3.6.1 # sybilSBML is necessary
# * DONE (CHNOSZ)
# * DONE (jsonlite)
module load cplex/12.9.0-academic  


cd /srv/scratch/z5434493/Am_queenslandica


# New:
/srv/scratch/z5434493/Softwares/gapseq/gapseq_20220825/gapseq find -p all -t Bacteria -b 200 -c 70 -l MetaCyc -y ./MAGs/AqS4.fasta

#√   -p keywords such as pathways or subsystems (for example amino,nucl,cofactor,carbo,polyamine)
#   -e Search by ec numbers (comma separated)
#   -r Search by enzyme name (colon separated)
#   -d Database: vmh or seed (default: seed)
#√   -t Taxonomic range for sequences to be downloaded (default: Bacteria)
#√   -b Bit score cutoff for local alignment (default: 200)
#   -i Identity cutoff for local alignment (default: 0)
#√   -c Coverage cutoff for local alignment (default: 75)
#   -s Strict candidate reaction handling (do _not_ use pathway completeness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: false)
#   -u Suffix used for output files (default: pathway keyword)
#   -a blast hits back against uniprot enzyme database
#√   -n Consider superpathways of metacyc database
#√   -l Select the pathway database (MetaCyc(2712), KEGG(523), SEED(666), all(3922); default: metacyc,custom)
#   -o Only list pathways found for keyword; default false)
#   -x Do not blast only list pathways, reactions and check for available sequences; default false
#   -q Include sequences of hits in log files; default false
#   -v Verbose level, 0 for nothing, 1 for pathway infos, 2 for full (default 1)
#   -k Do not use parallel
#   -g Exhaustive search, continue blast even when cutoff is reached (default false)
#   -z Quality of sequences for homology search: 1:only reviewed (swissprot), 2:unreviewed only if reviewed not available, 3:reviewed+unreviewed, 4:only unreviewed (default 2)
#   -m Limit pathways to taxonomic range (default )
#NEW -w Use additional sequences derived from gene names (default true)
#NEW -y Print annotation genome coverage (default false)
#NEW -j Quit if output files already exist (default false)
#NEW -U Do not use gapseq sequence archive and update sequences from uniprot manually (very slow) (default false)


# New:
/srv/scratch/z5434493/Softwares/gapseq/gapseq_20220825/gapseq find-transport -b 100 -c 70 ./MAGs/AqS4.fasta

#√   -b bit score cutoff for local alignment (default: 50)
#   -i identity cutoff for local alignment (default: 0)
#√   -c coverage cutoff for local alignment (default: 75)
#   -q Include sequences of hits in log files; default false
#   -k do not use parallel
#   -m only check for this keyword/metabolite (default )

# # # -b neg
# New:
/srv/scratch/z5434493/Softwares/gapseq/gapseq_20220825/gapseq draft -r ./AqS4-all-Reactions.tbl -t ./AqS4-Transporter.tbl -p ./AqS4-all-Pathways.tbl -c ./AqS4.fasta -u 100 -l 50 -b neg
# a few reactions should be added/modified to created draft model!!!

#√     -r|--blast.res          Blast-results table generated by gapseq.sh.
#√     -t|--transporter.res    Blast-results table generated by transporter.sh.
#√     -b|--biomass            Gram "pos" OR "neg" OR "archae" OR "auto"? Default: "auto". Please note: if set to "auto", the external programms barrnap, usearch, and bedtools are required.
#     -n|--model.name         Name of draft model network. Default: the basename of "blast.res"
#     -c|--genome.seq         If gram is set to "auto", the genome sequence is required to search for 16S genes, which are used to predict gram-staining.
#√     -u|--high.evi.rxn.BS    Reactions with an associated blast-hit with a bitscore above this value will be added to the draft model as core reactions (i.e. high-sequence-evidence reactions)
#√     -l|--min.bs.for.core    Reactions with an associated blast-hit with a bitscore below this value will be considered just as reactions that have no blast hit.
#     -o|--output.dir         Directory to store results. Default: "." (alternatives not yet implemented)
#     -s|--sbml.output        Should the gapfilled model be saved as sbml? Default: FALSE (export not yet implemented)
#     -p|--pathway.pred       Pathway-results table generated by gapseq.sh.
#     -a|--curve.alpha        Exponent coefficient for transformation of bitscores to reaction weights for gapfilling. (Default: 1 (neg-linear))

# Medium prediction
/srv/scratch/z5434493/Softwares/gapseq/gapseq_20220825/gapseq medium -m AqS4-draft.RDS -p AqS4-all-Pathways.tbl


# module unload R/3.6.1 # sybilSBML is necessary
# module add R/3.5.3 # sybilSBML is necessary
# 
/srv/scratch/z5434493/Softwares/gapseq/gapseq_20220825/gapseq fill -m ./AqS4-draft.RDS -n AqS4-medium.csv -c ./AqS4-rxnWeights.RDS -g ./AqS4-rxnXgenes.RDS -b 50 -o ./2022_pred_diet_AqS4 -r TRUE
# New:
# 1. sub2pwy.csv add 4 compounds: Oxidized-ferredoxins, Reduced-ferredoxins, Oxidized-Plastocyanins and Plastocyanin-Reduced
# 2. diet with single carbon source to seawater?

#√     -m|--model                  GapSeq-Draft-Model to be gapfilled (RDS or SBML)
#     -h|--help                   help
#√     -n|--media                  tab- or komma separated table for media components. Requires three named columns: 1 - "compounds" (for metab. IDs), 2 - "name" (metab. name), 3 - "maxFlux" (maximum inflow flux)
#     -f|--full.model             RDS file of the full (dummy) model. (ask Silvio for it :) ). Defaut: dat/full.model.RDS
#     -t|--target.metabolite      ID (without compartment suffix) of metabolite that shall be produced. Default: cpd11416 (Biomass)
#√     -c|--rxn.weights.file       Reaction weights table generated by gapseq function "generate_GSdraft.R" (RDS format).
#√     -g|--rxnXgene.table         Table with gene-X-reaction associations as generated by the "generate_GSdraft.R" (RDS format)
#√     -b|--bcore                  Minimum bitscore for reaction associated blast hits to consider reactions as core/candidate reactions. Default: 50
#√     -o|--output.dir             Directory to store results. Default: "gapfill".
#     -s|--sbml.output            Should the gapfilled model be saved as sbml? Default: FALSE
#     -q|--quick.gf               perform only step 1 and 2. Default: FALSE
#     -l|--limit                  Test metabolite to which search is limitted
#     -x|--no.core                Use always all reactions insteadof core reactions with have sequence evidence. Default: FALSE
#     -v|--verbose                Verbose output and printing of debug messages. Default: FALSE
#√     -r|--relaxed.constraints    Save final model as unconstraint network (i.e. all exchange reactions are open). Default: FALSE
