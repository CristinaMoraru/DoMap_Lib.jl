
export set_ProjReadmapIndexing, set_ProjReadMapMapping, ProjDoMap_fun

const ALLOWED_DoMAP = Dict(
    "readtypes" => ("short", "long", "mixed"),
    "maptool" => ("bowtie2", "minimap2"),
    "tokeep" => ("sam", "bam", "sam-bam")
)

function set_ProjReadmapIndexing(pd::String, args::Vector{String}, maptool::String, inref::FnaP, num_threads::Int64)
    if maptool == "bowtie2"
        #region bowtie indexing
            bowtie2index_args = [
                "projtype=index",
                "pd=$pd",         # here the pd is the same as the one for the whole proj, because DoBowtie is a complex workflow, with multiple steps
                "inref=$(inref.p)",
                "cpu=$num_threads",
                "rm_prev=false"
                ]
            index_proj = [initialize_bowtieproj(bowtie2index_args)]
        #endregion
    elseif maptool == "minimap2"  # if the assembly is mixed, then I need two indexes, one for short reads and one for long reads
        minimap2_p = extract_inPaths(args, "minimap2_p")

        #region minimap2 indexing
        if readtypes in ["short", "mixed"]
            MiniMap2Index_Ds = "$(pd)/minimap2index/short"
            MiniMap2IndexLogs_Ds = "$(MiniMap2Index_Ds)/logs"
            #MiniMap2IndexLogs_Ds = "$(minimap2Logs_D)/index/short"
            my_mkpath([MiniMap2Index_Ds, MiniMap2IndexLogs_Ds])
            index_projs = ProjMiniMap2Index(MiniMap2Index_Ds, WrapCmd(; cmd = RunMiniMap2IndexCmd(minimap2_p, "sr", inref, MmiP("$(MiniMap2Index_Ds)/$(indexname).mmi")), 
                log_p = "$(MiniMap2IndexLogs_Ds)/log_index.txt", err_p = "$(MiniMap2IndexLogs_Ds)/err_index.txt", exit_p = "$(MiniMap2IndexLogs_Ds)/exit_index.txt"))#, env = minimap2_env))
        end

        if readtypes in ["long", "mixed"]
            MiniMap2Index_Dl = "$(pd)/minimap2index/long"
            MiniMap2IndexLogs_Dl = "$(MiniMap2Index_Dl)/logs"
            my_mkpath([MiniMap2Index_Dl, MiniMap2IndexLogs_Dl])
            index_projl = ProjMiniMap2Index(MiniMap2Index_Dl, WrapCmd(; cmd = RunMiniMap2IndexCmd(minimap2_p, "map-ont", inref, MmiP("$(MiniMap2Index_Dl)/$(indexname).mmi")), 
                log_p = "$(MiniMap2IndexLogs_Dl)/log_index.txt", err_p = "$(MiniMap2IndexLogs_Dl)/err_index.txt", exit_p = "$(MiniMap2IndexLogs_Dl)/exit_index.txt"))#, env = minimap2_env))
        end

        if readtypes == "short"
            index_proj = [index_projs]
        elseif readtypes == "long"
            index_proj = [index_projl]
        elseif readtypes == "mixed"
            index_proj = [index_projs, index_projl]
        end
        #endregion
    end

    return index_proj
end

function set_ProjReadMapMapping(pd::String, args::Vector{String}, maptool::String, inref::FnaP, num_threads::Int64, index_proj::Vector{T}, inreadpairs::DataFrame) where T <: ProjReadMapIndexing
    samtools_p = extract_inPaths(args, "samtools_p")

    tokeep = extract_args(args, "tokeep", "bam"; allowed = ALLOWED_DoMAP["tokeep"])

    readpair_vec = Vector{Union{Missing, PerReadPair}}(missing, nrow(inreadpairs))

    samview_flag4exclusion = extract_args(args, "samview_flag4exclusion", Bool, "false")
    if samview_flag4exclusion == true
        samview_flag4exclusion_val = extract_args(args, "samview_flag4exclusion_val")
    else
        samview_flag4exclusion_val = ""
    end

    for row in 1:nrow(inreadpairs)

        read1 = inreadpairs[row, 1]
        read2 = inreadpairs[row, 2]
        if read2 == ""
            readname= read1
        else
            readname = stringIntersect(read1, read2)[1]
            # remove the "interleave...." from readname, leave only the sample name
        end

        if maptool == "bowtie2"
            #region bowtie2 mapping
            #shrinksam_p = extract_inPaths(args, "shrinksam_p")

            # bowtie2 mapping obj for each row in the Df
            bowtie2map_args = [
                "projtype=map",
                "pd=$pd/map/$readname", 
                "indexD=$(index_proj[1].BowtieBuild.cmd.indexD)",
                "indexname=$(index_proj[1].BowtieBuild.cmd.indexname)",
                "inref=$(inref.p)",
                "read1=$(read1)",
                "read2=$(read2)",
                "preset_endtoend=sensitive",
                "cpu=$num_threads",
                #"shrinksam_p=$shrinksam_p",
                "samtools_p=$samtools_p",
                "rm_prev=false",
                "samview_flag4exclusion=$samview_flag4exclusion",
                "samview_flag4exclusion_val=$samview_flag4exclusion_val"
                ]

            map_proj = initialize_bowtieproj(bowtie2map_args)
        #endregion
        elseif maptool == "minimap2"
            minimap2_p = extract_inPaths(args, "minimap2_p")

            #region minimap2 mapping
            minimap_parm = inreadpairs[row, 3]
            MiniMap2Alig_D = "$(pd)/map/minimap2/$readname"
            MiniMap2AligLogs_D = "$(MiniMap2Alig_D)/logs"
            my_mkpath([MiniMap2Alig_D, MiniMap2AligLogs_D])

            outbase = "$MiniMap2Alig_D/$(indexname)_$(readname)"
            out_f = "$(outbase).sam" |> SamP
            #outs_f = "$(outbase)_shrink.sam" |> SamP
            outview2sam_f = "$(outbase)_filt.sam" |> SamP
            outsorts_f = "$(outbase)_filt_sorted.sam" |> SamP
            outsortb_f = "$(outbase)_filt_sorted.bam" |> BamP
            covout_p = "$(outbase)_filt_sorted_coverage.tsv" |> TableP

            # minimap2
            if minimap_parm == "sr"    
                MiniMap2Alig = WrapCmd(cmd = RunMiniMap2AligCmd(minimap2_p, num_threads, "sr", index_projs.MiniMap2Index.cmd.index, FastaQP(read1), FastaQP(read2), out_f),
                    log_p = "$(MiniMap2AligLogs_D)/log.txt", err_p = "$(MiniMap2AligLogs_D)/err.txt", exit_p = "$(MiniMap2AligLogs_D)/exit.txt")#, env = minimap2_env)
            elseif minimap_parm == "map-ont"
                MiniMap2Alig = WrapCmd(cmd = RunMiniMap2AligCmd(minimap2_p, num_threads, "map-ont", index_projl.MiniMap2Index.cmd.index, FastaQP(read1), missing, out_f),
                    log_p = "$(MiniMap2AligLogs_D)/log.txt", err_p = "$(MiniMap2AligLogs_D)/err.txt", exit_p = "$(MiniMap2AligLogs_D)/exit.txt")#, env = minimap2_env)
            end
            
            
            # removes unpaired reads
            #=Shrinksam = WrapCmd(; cmd = RunShrinksamCmd(shrinksam_p, out_f, outs_f, false), 
                log_p = "$MiniMap2AligLogs_D/logShrink.txt", err_p = "$MiniMap2AligLogs_D/errShrink.txt", exit_p = "$MiniMap2AligLogs_D/exitShrink.txt")=#

            # removes reads with Flag 3584 (read fails quality check, read is PCR or optical duplicate, supplementary alignment; keeps the header)
            if samview_flag4exclusion == false
                more_opts = ["-h"]
            else
                more_opts = ["-h", "-F", "$samview_flag4exclusion_val"]
            end

            SamtoolsView2Sam = WrapCmd(; cmd = RunSamtoolsViewCmd(samtools_p, more_opts, num_threads, out_f, outview2sam_f),              # now the input for this step is the original sam file, not the shrinked one
                log_p = "$MiniMap2AligLogs_D/logSamtoolsView2Sam.txt", err_p = "$MiniMap2AligLogs_D/errSamtoolsView2Sam.txt", exit_p = "$MiniMap2AligLogs_D/exitSamtoolsView2Sam.txt")
                
            # sorts the sam file
            SamtoolsSort2Sam = WrapCmd(; cmd = RunSamtoolsSortCmd(samtools_p, ["-O", "SAM"], num_threads, outview2sam_f, outsorts_f),
                log_p = "$MiniMap2AligLogs_D/logSamtoolsSort2Sam.txt", err_p = "$MiniMap2AligLogs_D/errSamtoolsSort2Sam.txt", exit_p = "$MiniMap2AligLogs_D/exitSamtoolsSort2Sam.txt")

            # converts the ordered sam to bam (keeps the header)
            SamtoolsView2Bam = WrapCmd(; cmd = RunSamtoolsViewCmd(samtools_p, ["-h", "-b"], num_threads, outsorts_f, outsortb_f),
                log_p = "$MiniMap2AligLogs_D/logSamtoolsView2Bam.txt", err_p = "$MiniMap2AligLogs_D/errSamtoolsView2Bam.txt", exit_p = "$MiniMap2AligLogs_D/exitSamtoolsView2Bam.txt")

            # calculates the coverage (ignores reads smaller than 50 bases)
            SamtoolsCoverage = WrapCmd(; cmd = RunSamtoolsCovCmd(samtools_p, ["-l", "50"], outsorts_f, covout_p),  #["--min-MQ", "30"] removes reads with low MAPQ score
                log_p = "$(MiniMap2AligLogs_D)/logSamtoolsCov.txt", err_p = "$(MiniMap2AligLogs_D)/errSamtoolsCov.txt", exit_p = "$(MiniMap2AligLogs_D)/exitSamtoolsCov.txt")

            map_proj = ProjMiniMap2Alig(MiniMap2Alig_D, MiniMap2Alig, SamtoolsView2Sam, SamtoolsSort2Sam, SamtoolsView2Bam, SamtoolsCoverage) #Shrinksam, 
            #endregion
        end

        # cov tables for MaxBin2 and Metabat2
        cov4MaxBin2_p = "$(map_proj.pd)/maxbin2_cov.tsv" |> TableP

        readpair_vec[row] = PerReadPair(readname, map_proj, cov4MaxBin2_p, tokeep)
    end
    
    return readpair_vec
end        

function ProjDoMap_fun(args::Vector{String})
    # num_threads
    num_threads = extract_args(args, "num_threads", Int64, 1, 1, 40)

    # continue or not
    cont = extract_args(args, "continue", Bool, "false")

    # inref
    inref = extract_inFiles(args, "inref", BioS_Gen.ALLOWED_EXT["FnaP"]) |> FnaP
    sampleName = getFileName(inref.p) 

    # folders
    pd_prefix = extract_args(args, "pd_prefix")
    pd = "$(pd_prefix)$(sampleName)"

    if cont == false
        do_pd(pd)
        dosteps = Dict(
            "indexing" => WorkflowStatus("do", "not_done"),
            "mapping" => WorkflowStatus("do", "not_done"))

        # reads
        inreadpairs_p = extract_inFiles(args, "inreadpairs", BioS_Gen.ALLOWED_EXT["TableP"]) |> TableP
        inreadpairs = CSV.read(inreadpairs_p.p, DataFrame; delim='\t', header=1)
        
        # params for readtypes and mapping program to use
        readtypes = extract_args(args, "readtypes", ALLOWED_DoMAP["readtypes"])
        maptool = extract_args(args, "maptool", ALLOWED_DoMAP["maptool"])

        if readtypes in ["long", "mixed"]
            maptool = "minimap2"           
        end

        proj = ProjSDoMap(pd = pd, sampleName = sampleName, cont = cont, maptool = maptool, readtypes = readtypes, inref = inref, inreadpairs = inreadpairs, dosteps = dosteps)
    else   # this section is in work, just a placehodler for the moment
        proj = load_proj("$(pd)/sproj.binary", pd, inref, sampleName)

        dosteps = proj.dosteps
    end

    mapindex = initialize_step(dosteps, "indexing", set_ProjReadmapIndexing, (pd, args, maptool, inref, num_threads), proj.mapindex, cont) 
    perreadpairs = initialize_step(dosteps, "mapping", set_ProjReadMapMapping, (pd, args, maptool, inref, num_threads, mapindex, inreadpairs), proj.perreadpairs, cont)

    proj = nothing
    sproj = ProjSDoMap("singleworkflow", pd, sampleName, cont, maptool, readtypes, inref, inreadpairs, dosteps, mapindex, perreadpairs)

    serialize("$(pd)/sproj.binary", sproj)
    return sproj
end