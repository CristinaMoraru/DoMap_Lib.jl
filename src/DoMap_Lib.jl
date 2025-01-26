module DoMap_Lib

using Revise
using CSV
using DataFrames
using Serialization

using BioS_Gen
using BioS_ProjsWFs
import BioS_ProjsWFs.run_workflow
import BioS_ProjsWFs.run_workflow!
using BioS_ExtTools


include("DoMap_Lib_00_DT_structs.jl")
include("DoMap_Lib_00_DT_MainProjConstructor.jl")
include("DoMap_Lib_000_Run.jl")


#args = ARGS

#=ProjMultiWorkflow - a multiple binning workflow
const args = [
    "projtype=multipleworkflow",
    "spd=/mnt/cephfs1/projects/SulfitobacterM53_ICBM16-18_transcriptome/test/outDoMap_4readpairs_ICBM16",
    "allrefs_params=/mnt/cephfs1/projects/SulfitobacterM53_ICBM16-18_transcriptome/01_Read_Mapping_DoMap/00_params/all_refs_params.tsv",
    "continue=false",
] =#


#=
const args = [
    "projtype=singleworkflow",
    "pd_prefix=/mnt/cephfs1/projects/SulfitobacterM53_ICBM16-18_transcriptome/test/outDoMap_2readpairs_2", #/data3/CLM_projs/RESIST_proposal/C13/mapping/out",     
    "inref=/data3/CLM_projs/SulfivirusesICBM_all/Sulfiviruses_72Fullgenomes/Sulfiviruses_72Full_individual/ICBM16.fasta",
    "inreadpairs=/mnt/cephfs1/projects/SulfitobacterM53_ICBM16-18_transcriptome/01_Read_Mapping_DoMap/00_params/read_pairs_ICBM16.tsv", #/data3/CLM_projs/RESIST_proposal/C13/mapping/params/reads_C13.tsv", # one row contains 
    
    "readtypes=short",
    "maptool=bbmap", #or bowtie2, minimap2 or bbmap. It defaults to minimap2 if readtypes is "long" or "mixed"
    "tokeep=bam",
    "samview_flag4exclusion=false",
    "samview_flag4exclusion_val=3584",

    "bbmap_p=/software/conda/soft/BBTools/bbmap_39.14/bbmap.sh",
    "minimap2_p=/software/conda/soft/minimap2-2.26_x64-linux/minimap2",
    "shrinksam_p=/software/conda/soft/shrinksam/shrinksam",
    "samtools_p=/software/conda/soft/samtools/samtools-1.18Inst/bin/samtools",

    #postmap filtering bbmap
    "idfilter=1.0",
    "subfilter=0",
    "insfilter=0",
    "delfilter=0",
    "indelfilter=0",
    "inslenfilter=0",
    "dellenfilter=0",
    "nfilter=0",
    "ambiguous=toss", 

    "continue=false",
    "num_threads=20"
] =#

#= To inactivate when used just as library
if "--help" in args
    println("DoMap - a workflow for read mapping with a choice of three aligners: BBMap, Bowtie or MiniMap2.")
end

println("Start DoMap!")

proj = initialize_workflow(args, ProjSDoMap_fun) # args


if proj.projtype == "singleworkflow"
    run_workflow(proj)
elseif proj.projtype == "multipleworkflow"
    run_workflowDoMap(proj)
end

println("DoMap is done!")
=#

end # module DoMap_Lib

