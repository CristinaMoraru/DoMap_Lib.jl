module prep_reads
using CSV
using DataFrames

read1 = Vector{String}()
read2 = Vector{String}()

for i in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "21", "22", "23", "24", "29", "30", "31", "32", "37", "38", "45", "46", "47", "53", "54", "61", "62"]
    push!(read1, "/data3/CLM_projs/RESIST_proposal/ExStreamIMP-31metagenomes/raw/ES22_IMP_S_C$(i)_interleaved_trim_clean_no_dedup.PE.1.fastq.gz")
    push!(read2, "/data3/CLM_projs/RESIST_proposal/ExStreamIMP-31metagenomes/raw/ES22_IMP_S_C$(i)_interleaved_trim_clean_no_dedup.PE.2.fastq.gz")
end

df = DataFrame(read1 = read1, read2 = read2)

out_p = "/data3/CLM_projs/RESIST_proposal/ExStreamIMP-31metagenomes/mapping/params/ExStreamIMP_31_in_reads.tsv" 
CSV.write(out_p, df, delim='\t', header=true)

end #module