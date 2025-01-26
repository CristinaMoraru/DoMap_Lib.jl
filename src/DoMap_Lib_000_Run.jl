export run_workflow_postmap, run_workflow_MiniMap2Alig, run_workflow_BBmapIndex, run_workflow_BBMapMap, run_workflow, run_workflowDoMap


function run_workflow_postmap(proj::Union{ProjMiniMap2Alig, ProjBBmapMap}, parentD::String)
    println("Starting post-mapping operations: filtering, sorting, converting to bam and coverage calculations.")

    #do_cmd(proj.Shrinksam, "Shrinksam", false)
    do_cmd(proj.SamtoolsView2Sam, "SamtoolsView2Sam", false, parentD)
    do_cmd(proj.SamtoolsSort2Sam, "SamtoolsSort", false, parentD)
    do_cmd(proj.SamtoolsView2Bam, "SamtoolsView", false, parentD)
    do_cmd(proj.SamtoolsCoverage, "SamtoolsCoverage", false, parentD)
end

function run_workflow_MiniMap2Alig(proj::ProjMiniMap2Alig, parentD::String)
    # map
    println("Starting minimap2 alignment.")
    do_cmd(proj.MiniMap2Alig, "minimap2", false, parentD)

    # postprocessing
    run_workflow_postmap(proj, parentD)

    return nothing
end

function run_workflow_BBmapIndex(proj::ProjBBMapIndex, parentD::String)
    #println("Starting BBMap Indexing")
    cd("$parentD/$(proj.pd)")
    
    do_cmd(proj.BBmapIndex, "BBMap indexing", false, parentD)

    return nothing
end

function run_workflow_BBMapMap(proj::ProjBBmapMap, parentD::String)
    # mapping
    #println("Starting BBMap mapping")
    cd("$(parentD)/$(proj.pd)")
    do_cmd(proj.BBMapMap, "BBMap mapping", false, parentD)

    # postprocessing
    run_workflow_postmap(proj, parentD)

    return(nothing)
end

function run_workflow(proj::ProjSDoMap)

    # indexing
    if proj.dosteps["indexing"].signal == "do"
        set2running!("indexing", proj)

        if proj.maptool == "bowtie2"
            # bowtie2 indexing
            run_workflow(proj.mapindex[1])  # here the index will always be 1, bc bowtie2 produces only one index
        elseif proj.maptool == "minimap2"
            # minimap2 indexing
            for i in eachindex(proj.mapindex)
                do_cmd(proj.mapindex[i].MiniMap2Index, "minimap2 indexing", false, proj.pd)   # here the index will be i, bc minimap2 can produce three types (each in one folder) of indexes
            end
        elseif proj.maptool == "bbmap"
            # bbmap indexing
            for i in eachindex(proj.mapindex)
                run_workflow_BBmapIndex(proj.mapindex[1], proj.pd) # here the index will always be 1, bc bbmap produces only one index
            end
        end
        
        set2finished!("indexing", proj)
    end


    # mapping
    if proj.dosteps["mapping"].signal == "do"
        set2running!("mapping", proj)

        for row in eachindex(proj.perreadpairs)
            if proj.maptool == "bowtie2"
                # bowtie2 mapping and coverage calculation
                run_workflow(proj.perreadpairs[row].mapping)
            elseif proj.maptool == "minimap2"
                # minimap2 mapping and coverage calculation
                run_workflow_MiniMap2Alig(proj.perreadpairs[row].mapping, proj.pd)
            elseif proj.maptool == "bbmap"
                # bbmap mapping and feature counts
                run_workflow_BBMapMap(proj.perreadpairs[row].mapping, proj.pd)
            end

            # generate abundace tables for MaxBin2 from SamtoolsCoverage data
            #genAbundMaxBin(proj.perreadpairs[row])

            # remove first two sam files
            println("Removing .sam files.")
            if proj.maptool == "bowtie2"
                rm_path([#proj.perreadpairs[row].mapping.Shrinksam.cmd.output.p,
                        "$(proj.pd)/$(proj.perreadpairs[row].mapping.SamtoolsView2Sam.cmd.out_f.p)"])
            elseif proj.maptool == "minimap2"
                rm_path(["$(proj.pd)/$(proj.perreadpairs[row].mapping.MiniMap2Alig.cmd.out_p.p)",
                        "$(proj.pd)/$(proj.perreadpairs[row].mapping.SamtoolsView2Sam.cmd.out_p.p)"])
            elseif proj.maptool == "bbmap"
                rm_path(["$(proj.pd)/$(proj.perreadpairs[row].mapping.BBMapMap.cmd.out_f.p)",
                        "$(proj.pd)/$(proj.perreadpairs[row].mapping.SamtoolsView2Sam.cmd.out_f.p)"])
            end

            #remove the last sam file if the user wants to keep only bam files
            if proj.perreadpairs[row].tokeep == "bam"
                rm_path(["$(proj.pd)/$(proj.perreadpairs[row].mapping.SamtoolsSort2Sam.cmd.out_f.p)"])
            end
        end

        # remove the index folder
        if proj.maptool == "bowtie2"
            rm_path(["$(proj.pd)/$(proj.mapindex[1].BowtieBuild.cmd.indexD)"])   #here the index is always 1, bc bowtie2 produces only one index
        elseif proj.maptool == "minimap2"
            for i in eachindex(proj.mapindex)
                rm_path(["$(proj.pd)/$(proj.mapindex[i].pd)"])
            end
        elseif proj.maptool == "bbmap"
            rm_path(["$(proj.pd)/$(proj.mapindex[1].pd)/ref"])
        end

        set2finished!("mapping", proj)
    end

end


function run_workflowDoMap(mproj::ProjMultiWorkflow)

    for i in eachindex(mproj.allSingleWorkflows)
        if mproj.dosteps[i].progress in ["not_done", "running"]
            println("
            ---------------------- Starting TO RUN the DoMap workflow for input file $i -------------------------
            ")
            
            set2running!(i, mproj)

            run_workflow(mproj.allSingleWorkflows[i])
           
            set2finished!(i, mproj)

            println("
            ---------------------- Finished RUNNING the DoMap workflow for input file $i ------------------------
            ")
        else
            println("
            ---------------------- Skipping TO RUN the DoMap workflow for input file $i, because its status in a previous run was $(mproj.dosteps[i].progress) ------------------------
            ")
        end
    end

    return nothing
end