module DpMap_pre_params
using CSV
using DataFrames

df_p = "/data3/CLM_projs/RESIST_proposal/ExStreamIMP-31metagenomes/mapping/params/ExStreamIMP_31_inref_paths.tsv" #"/data3/CLM_projs/AcBaMe/in/acbame_contigs.tsv" #"/data3/CLM_projs/Nitrospira_proj/MeatAOB_NOB_Seb_Luck_analApril2024/viruspred/in/params/sebluck_12met_dovip_in.tsv"
df = CSV.read(df_p, DataFrame; delim='\t', header=true)

df[!, :projtype] = fill("singleworkflow", nrow(df))
df[!, :continue] = fill("false", nrow(df))


df[!, :inreadpairs] = fill("/data3/CLM_projs/RESIST_proposal/ExStreamIMP-31metagenomes/mapping/params/ExStreamIMP_31_reads.tsv", nrow(df))
df[!, :readtypes] = fill("short", nrow(df))
df[!, :maptool] = fill("bowtie2", nrow(df))
df[!, :tokeep] = fill("bam", nrow(df))

df[!, :minimap2_p] = fill("/software/conda/soft/minimap2-2.26_x64-linux/minimap2", nrow(df))
df[!, :shrinksam_p] = fill("/software/conda/soft/shrinksam/shrinksam", nrow(df))
df[!, :samtools_p] = fill("/software/conda/soft/samtools/samtools-1.18Inst/bin/samtools", nrow(df))

df[!, :samview_flag4exclusion] = fill("fase", nrow(df))

#df[!, :] = fill(, nrow(df))
#df[!, :] = fill(, nrow(df))

df[!, :num_threads] = fill("15", nrow(df))

out_p = "/data3/CLM_projs/RESIST_proposal/ExStreamIMP-31metagenomes/mapping/params/ExStreamIMP_31_in_params.tsv" #/data3/CLM_projs/AcBaMe/in/acbame_contigs_dovip_params.tsv"
CSV.write(out_p, df, delim='\t', header=true)


end # module