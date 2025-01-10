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

#ProjMultiWorkflow - a multiple binning workflow
const args = [
    "projtype=multipleworkflow",
    "spd=/data3/CLM_projs/RESIST_proposal/C13/mapping/out/all",
    "allrefs_params=/data3/CLM_projs/RESIST_proposal/C13/mapping/params/all_params.tsv",
    "continue=false",
] #


#=
const args = [
    "projtype=singleworkflow",
    "pd_prefix=/data3/CLM_projs/TEST_Workflows/outDoMap_", #/data3/CLM_projs/RESIST_proposal/C13/mapping/out",
    "inref=/data3/CLM_projs/TEST_Workflows/inDoMap/AOB2Aal2016_NonIntegrated_outer_thresholding.fna", #/data3/CLM_projs/RESIST_proposal/C13/viruses/ExStreamIMP_C13_ALL_NonIntegrated_outer_thresholding.fna",
    "inreadpairs=/data3/CLM_projs/TEST_Workflows/inDoBiP/params/readpairs_2ref.tsv", #/data3/CLM_projs/RESIST_proposal/C13/mapping/params/reads_C13.tsv", # one row contains 
    
    "readtypes=short",
    "maptool=bowtie2", #or minimap2' defaults to minimap2 if readtypes is "long" or "mixed"
    "tokeep=bam",
    "samview_flag4exclusion=false",
    "samview_flag4exclusion_val=3584",

    "minimap2_p=/software/conda/soft/minimap2-2.26_x64-linux/minimap2",
    "shrinksam_p=/software/conda/soft/shrinksam/shrinksam",
    "samtools_p=/software/conda/soft/samtools/samtools-1.18Inst/bin/samtools",

    "continue=false",
    "num_threads=30"
] =#

if "--help" in args
    println("DoViP - a workflow for virus prediction in metagenomes.")
end

println("Start DoMap!")

proj = initialize_workflow(args, ProjDoMap_fun) # args
Threads.nthreads() = 20

if proj.projtype == "singleworkflow"
    run_workflow(proj)
elseif proj.projtype == "multipleworkflow"
    run_workflowDoMap(proj)
end

println("DoMap is done!")

end # module DoMap_Lib


#
#= ProjMultiWorkflow - a multiple binning workflow
const args = [
    "projtype=multipleworkflow",
    "spd=/mnt/XIO_2/data3/CLM_projs/TEST_Workflows/outDoBiP1",
    "allrefs_params=/mnt/XIO_2/data3/CLM_projs/TEST_Workflows/inDoBiP/params/master_2refs.tsv",
    "continue=false",
]=#