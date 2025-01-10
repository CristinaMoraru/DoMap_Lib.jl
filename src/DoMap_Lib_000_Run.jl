export run_workflow

function run_workflow(proj::ProjMiniMap2Alig)
    # map
    println("Starting minimap2 alignment.")
    do_cmd(proj.MiniMap2Alig, "minimap2", false)

    # postprocessing
    println("Starting post-mapping operations: shrinking, filtering, sorting and converting to bam.")
    #do_cmd(proj.Shrinksam, "Shrinksam", false)
    do_cmd(proj.SamtoolsView2Sam, "SamtoolsView2Sam", false)
    do_cmd(proj.SamtoolsSort2Sam, "SamtoolsSort", false)
    do_cmd(proj.SamtoolsView2Bam, "SamtoolsView", false)
    do_cmd(proj.SamtoolsCoverage, "SamtoolsCoverage", false)
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
                do_cmd(proj.mapindex[i].MiniMap2Index, "minimap2 indexing", false)   # here the index will be i, bc minimap2 can produce three types (each in one folder) of indexes
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
                run_workflow(proj.perreadpairs[row].mapping)
            end

            # generate abundace tables for MaxBin2 from SamtoolsCoverage data
            #genAbundMaxBin(proj.perreadpairs[row])

            # remove first two sam files
            println("Removing .sam files.")
            if proj.maptool == "bowtie2"
                rm_path([#proj.perreadpairs[row].mapping.Shrinksam.cmd.output.p,
                        proj.perreadpairs[row].mapping.SamtoolsView2Sam.cmd.out_f.p])
            else proj.maptool == "minimap2"
                rm_path([proj.perreadpairs[row].mapping.MiniMap2Alig.cmd.out_p.p,
                        proj.perreadpairs[row].mapping.SamtoolsView2Sam.cmd.out_f.p])
            end

            #remove the last sam file if the user wants to keep only bam files
            if proj.perreadpairs[row].tokeep == "bam"
                rm_path([proj.perreadpairs[row].mapping.SamtoolsSort2Sam.cmd.out_f.p])
            end
        end

        # remove the index folder
        if proj.maptool == "bowtie2"
            #for i in eachindex(proj.mapindex)
                rm_path([proj.mapindex[1].BowtieBuild.cmd.indexD])   #here the index is always 1, bc bowtie2 produces only one index
            #end
        elseif proj.maptool == "minimap2"
            for i in eachindex(proj.mapindex)
                rm_path([proj.mapindex[i].pd])
            end
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