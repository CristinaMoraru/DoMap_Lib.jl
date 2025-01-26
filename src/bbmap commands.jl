[bbmap_p="/software/conda/soft/BBTools/bbmap_39.14/bbmap.sh"]

#= cd work directory to the folder for the reference?
     - the   =# 

#write the index, only once per workflow, in the folder of the reference named /ref/
bbmap.sh ref=A.fa   # paste0(BBMAP_exe, " ref=", genomes_path_ls, " -Xmx4g k=13")

#when doing the mapping, make sure to cd to the folder of the sample, in which the ref subfolder exists and contains the reference
# save bbmap data in bbmap_out
bbmap.sh in=reads.fq