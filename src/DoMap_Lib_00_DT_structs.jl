#export ProjReadMap, 
export ProjMiniMap2Index, ProjMiniMap2Alig, ProjBBMapIndex, ProjBBmapMap, PerReadPair, ProjSDoMap

#region Minimap
struct ProjMiniMap2Index <: ProjReadMapIndexing
    pd::String
    MiniMap2Index::WrapCmd{RunMiniMap2IndexCmd}
end

struct ProjMiniMap2Alig <: ProjReadMapMapping
    pd::String
    MiniMap2Alig::WrapCmd{RunMiniMap2AligCmd}
    #Shrinksam::WrapCmd{RunShrinksamCmd}
    SamtoolsView2Sam::WrapCmd{RunSamtoolsViewCmd}
    SamtoolsSort2Sam::WrapCmd{RunSamtoolsSortCmd}
    SamtoolsView2Bam::WrapCmd{RunSamtoolsViewCmd}
    SamtoolsCoverage::WrapCmd{RunSamtoolsCovCmd}
end
#endregion

#region BBMap
struct ProjBBMapIndex <: ProjReadMapIndexing
    pd::String
    BBmapIndex::WrapCmd{RunBBMapIndex}
end

struct ProjBBmapMap <: ProjReadMapMapping
    pd::String
    BBMapMap::WrapCmd{RunBBMapMap}
    SamtoolsView2Sam::WrapCmd{RunSamtoolsViewCmd}
    SamtoolsSort2Sam::WrapCmd{RunSamtoolsSortCmd}
    SamtoolsView2Bam::WrapCmd{RunSamtoolsViewCmd}
    SamtoolsCoverage::WrapCmd{RunSamtoolsCovCmd}
end

#endregion

#region MapProject
struct PerReadPair
    readname::String
    mapping::ProjReadMapMapping
    cov4MaxBin2::Union{Nothing, TableP}
    tokeep::String
end


Base.@kwdef mutable struct ProjSDoMap <: BioinfSProj
    projtype::String = "singleworkflow"
    pd::String
    sampleName::String
    cont::Bool
    maptool::String
    readtypes::String
    inref::FnaP
    inreadpairs::DataFrame
    dosteps::Dict{String, WorkflowStatus}
    mapindex::Union{Missing, Vector{T}} where T <: ProjReadMapIndexing = missing
    perreadpairs::Union{Missing, Vector{PerReadPair}} = missing
end
#endregion

