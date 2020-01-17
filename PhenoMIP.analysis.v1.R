# Version 6 of the MIP analysis will do a little more with the data in terms of how it spits things out
# Ideally there will be two types of data to analyse: mapping data, growth analysis data
# It will now produce 
# 1. an individual file for each experiment
# 2. a summary file of all the experiments
# 3. in anticipation of mapping, it will generate a map of the genome and outline the MIP ratio values per experiment
# 4. (140315) re-introduced analysis of the value files to specifically check the SNV at the given position. It will then remove samples that do not have
# a high enough value

# Needed updates:
# 140728: Should alter molecular barcode step to look for consensus on samples with repeat barcodes. If none, then remove both from the read set.
# This needs to be fixed in the get.MIP.barcodes() function
# 140819: v7 Iteration update. Includes looking at repeat barcodes (low percentage of them) and comparing the SNV between sets for conflict
# This update mainly affects function get.MIP.Bardcodes
# 150326: v9 New update started to re-evaluate samples and post-analysis and remove overly abundant strains from analysis (re-evaluate abundance of lower strains)
# 150406: v10 New update to deal with NextSeq data samples. When demultiplexed there are 4 lanes for each sample making the format [sample]_S[num]_L00[1-4]_R1_001.fastq
# 150914: v13 Update so that multiple SNV sites within a MIP can be evaluated and quality filtered.
# 160302: v15 Update to molecular barcode duplication (UMI). Based on random sequences, and from a lab meeting Chau just had, it may be prudent
# to count single-mismatch sequences for UMI as the same. Hiatt paper removes homopolymers >=4 in terms of UMI but doesn't mention any mismatches.
# 161209: v16 updated all grep uses to include perl=TRUE for faster speed
# 170408: v18 upgrade to use a make.fuzzy.query parameter for rDNA analysis. 

source("genomic.helper.R")

## ---------- Global variables ---------- ##

prog.ver = 18
MIP.read.length = 75 # This is the expected number of bases to be read into the program
MIP.barcode.length = 12 # This is the expected number of bases for the molecular barcode


#Chromosome information
LG.sizes = c(15072423, 15279345, 13783700, 17493793, 20924149, 17718866)

# Set the minimum quality score of a single basepair to 20? that's a 0.01% chance of being false. 30 = 0.001
min.quality = 30
# quality scores start at 14? /0123456789<>ABCDEFGH
output.info = "run.info.txt"
# maximum homopolymer length of unique molecular identifiers?
max.poly = 12

# Set the number of fuzzy bases allowed in our matching scheme
max.mismatch = 0


## ---------- Batch programs for running on data ---------- ##
# Use this if you have multiple sequencing data sets either MIP-MAP or PhenoMIP that need to be analysed at the same time
# See example batch file "batch.run.example.txt" for more information

batch.run = function(batch.file) {
  
  batch.data = read.table(batch.file, header=TRUE, sep="\t", colClasses = "character")  
  
  batch.results = list()
  
  for (i in 1:nrow(batch.data)) {
    
    if (batch.data$version[i] == "map") {
      
      batch.results[[i]] = batch.mapping(batch.data$fastq.file[i], batch.data$mapping.file[i], batch.data$output.dir[i], as.numeric(batch.data$start[i]), as.numeric(batch.data$end[i]), as.logical(batch.data$split[i]))
      
    } else if (batch.data$version[i] == "abundance") {
      
      batch.results[[i]] = batch.multiplex(batch.data$fastq.file[i], batch.data$mapping.file[i], batch.data$output.dir[i], as.numeric(batch.data$start[i]), as.numeric(batch.data$end[i]), as.logical(batch.data$split[i]))
      
    } 
        
  }
  
  return (batch.results)
  
}


# This is just a batch command to both split the fastq files and analyse their data for mapping
batch.mapping = function(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row, split, num.mismatch=0) {
  
  assign("max.mismatch", num.mismatch, envir=.GlobalEnv)
  
  #file.split.results = make.seq.files(MS.fastq.file, "seqFiles", "seqValues")
  if (isTRUE(split)) {
    make.seq.files(MS.fastq.file, "seqFiles", "seqValues", start.row, end.row)
  }
  
  mapping.results = analyse.mapping(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row)
  
  #return (list(file.split.results, mapping.results))
  return (mapping.results)
}

# This is the batch command to analyse multiple PhenoMIP data sets. 
# This has been replaced by the batch.run command
batch.batch.multiplex = function(batch.file) {
  
  batch.data = read.table(batch.file, header=TRUE, sep="\t", colClasses = "character")  
  
  batch.results = list()
  
  for (i in 1:nrow(batch.data)) {
    
    batch.results[[i]] = batch.multiplex(batch.data$fastq.file[i], batch.data$mapping.file[i], batch.data$output.dir[i], as.numeric(batch.data$start[i]), as.numeric(batch.data$end[i]))
    
  }
  
  return (batch.results)
}

# A batch command to look at multiplexed/pooled strain data
# Called within the batch.run function
batch.multiplex = function(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row, split) {
  
  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  assign("output.info", paste(c(output.directory, "run.info.txt"), collapse = ""), envir=.GlobalEnv)
  dir.create(output.directory)
  
  
  lapply(c(paste0("Running MIP.analysis.v",prog.ver, ".R at ", Sys.time())), write, output.info, append=FALSE, ncolumns=1000)
  lapply(c(paste0("min.quality: ", min.quality, "; max.poly: ", max.poly)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("max.mismatch: ", max.mismatch)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("MS.fastq.file: ", MS.fastq.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("mapping.info.file: ", mapping.info.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("output.directory: ", output.directory)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("start row: ", start.row)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("end row: ", end.row)), write, output.info, append=TRUE, ncolumns=1000)
  
  if (isTRUE(split)) {
    make.seq.files(MS.fastq.file, "seqFiles", "seqValues", start.row, end.row)
  }  
  
  multiplex.results = analyse.multiplex(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row)
  
  return (multiplex.results)
}

## ---------- Main programs for analysis ---------- ##

# This will take in fastq data from a mapping run and spit out
# 1. summary files for each experiment
# 2. line graph for each experiment 
analyse.mapping = function(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row) {
  
  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  #output.info = paste(c(output.directory, "run.info.txt"), collapse = "")
  assign("output.info", paste(c(output.directory, "run.info.txt"), collapse = ""), envir=.GlobalEnv)
  
  # Create the output directory
  dir.create(output.directory)
  
  # create a table of all the fastq file names
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")  
  # Trim down the table to what we're really interested in from it
  MS.fastq.data = MS.fastq.data[start.row:end.row,]
  
  # generate a list of all the fastq file locations
  MS.fastq.loc = MS.fastq.data$seqLocation
  
  # generate a list of all the fastq value locations
  MS.fastq.val = MS.fastq.data$valLocation
  
  # generate a list of output prefixes for the data
  output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(output.directory, x), collapse="")))
  
  # Read the mapping information into a table examination
  MIP.info = read.table(mapping.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  MIP.data = list()

  lapply(c(paste0("Running MIP.analysis.v",prog.ver, ".R at ", Sys.time())), write, output.info, append=FALSE, ncolumns=1000)
  lapply(c(paste0("min.quality: ", min.quality, "; max.poly: ", max.poly)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("max.mismatch: ", max.mismatch)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("MS.fastq.file: ", MS.fastq.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("mapping.info.file: ", mapping.info.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("output.directory: ", output.directory)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Table read complete at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Total experiments: ", (end.row-start.row+1))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Start row: ", (start.row))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("End row: ", (end.row))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Total MIPS to evaluate: ", nrow(MIP.info))), write, output.info, append=TRUE, ncolumns=1000)
  
  mapping.results = unlist(generate.graph.data(MIP.info$MIP.name, MIP.info$snv.chr, as.numeric(MIP.info$snv.loc)))
  read.summary = mapping.results
  # How many high-quality, gap-fill matches do we have?
  gf.match.summary = mapping.results
  fuzzy.match.summary = mapping.results
  
  # Start iterating through all the fastq files and analysing the MIP data
  for (i in 1:length(MS.fastq.loc)) {
  #for (i in start.row:end.row) {
    lapply(c(paste0("Starting analysis of ", MS.fastq.loc[i], " at ", Sys.time(), " with ", MS.fastq.data$seq.file.size[i], " reads.")), write, output.info, append=TRUE, ncolumns=1000)
    
    # Go grab the relevant information about how the MIPS did 
    MIP.data[[i]] = get.MIP.Barcodes(MIP.info, MS.fastq.loc[i], MS.fastq.val[i])    
    
    # Now cbind it to the mapping information to save later
    mapping.results = cbind(mapping.results, unlist(MIP.data[[i]]$snv.wt.unique.ratio))
    
    # With each iteration we should try to save data about each MIP for each experiment
    # For now the number of reads would be good to know
    read.summary = cbind(read.summary, unlist(MIP.data[[i]]$total.match.lig))
    gf.match.summary = cbind(gf.match.summary, (unlist(MIP.data[[i]]$wt.bc.unique)+ unlist(MIP.data[[i]]$snv.bc.unique)))
    fuzzy.match.summary = cbind(fuzzy.match.summary, (unlist(MIP.data[[i]]$wt.bc.fuzzy) + unlist(MIP.data[[i]]$snv.bc.fuzzy)))
    
    # Now save this table to the running directory      
    write.table(MIP.data[[i]], file=paste(c(output.prefix.list[i], "results", "txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)


    # 170222: Try to write more of the results as they are generated, in case the program crashes out due to large file sizes
  
    colnames(mapping.results) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    colnames(read.summary) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    colnames(gf.match.summary) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    colnames(fuzzy.match.summary) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    
    write.table(mapping.results, file=paste(c(output.directory, "mappingdata.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    write.table(read.summary, file=paste(c(output.directory, "read.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(gf.match.summary, file=paste(c(output.directory, "gap.fill.match.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    write.table(fuzzy.match.summary, file=paste(c(output.directory, "fuzzy.match.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    
    lapply(c(paste0("Analysis completed at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)  
  }
  
  colnames(mapping.results) = c("MIP.name", "location", MS.fastq.data$exp.name)
  colnames(read.summary) = c("MIP.name", "location", MS.fastq.data$exp.name)
  colnames(gf.match.summary) = c("MIP.name", "location", MS.fastq.data$exp.name)
  
  # Now collate all the data into a larger file of information
  all.MIP.info = summarize.MIPS(MIP.info, MIP.data)
  write.table(all.MIP.info, file=paste(c(output.directory, "MIP.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  write.table(mapping.results, file=paste(c(output.directory, "mappingdata.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  write.table(read.summary, file=paste(c(output.directory, "read.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  write.table(gf.match.summary, file=paste(c(output.directory, "gap.fill.match.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  
  normalized.mapping = normalize.data(paste(c(output.directory, "mappingdata.txt"), collapse=""), "./Normalization/150121.VC20019.norm.curv.txt", paste(c(output.directory, "norm.mappingdata.txt"), collapse=""))
  
  lapply(c(paste0("Analysis complete at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  return (list(MIP.data, all.MIP.info, mapping.results))
}

### ----------------------- analysis of multiplex population data for genetic interactions -----------------------###

# This will take in fastq data from a mapping run and spit out
# 1. summary files for each experiment
# 2. line graph for each experiment 
analyse.multiplex = function(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row) {
  
#   output.directory = paste(c("./", output.directory, "/"), collapse = "")
#   output.info = paste(c(output.directory, "run.info.txt"), collapse = "")
    
  # Create the output directory
  dir.create(output.directory)
  
  # create a table of all the fastq file names
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")  
  
  # generate a list of all the fastq file locations
  MS.fastq.loc = MS.fastq.data$seqLocation
  
  # generate a list of all the fastq value locations
  MS.fastq.val = MS.fastq.data$valLocation
  
  # generate a list of output prefixes for the data
  output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(output.directory, x), collapse="")))
  
  # Read the mapping information into a table examination
  MIP.info = read.table(mapping.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  MIP.data = list()

  lapply(c(paste0("Table read complete at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Total experiments: ", (end.row-start.row+1))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Start row: ", (start.row))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("End row: ", (end.row))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Total MIPS to evaluate: ", nrow(MIP.info))), write, output.info, append=TRUE, ncolumns=1000)
  
  mapping.results = cbind(MIP.info$MIP.name, MIP.info$snv.strain)
  
  read.summary = mapping.results
  
  total.rows = end.row - start.row + 1
  
  for (i in 1:total.rows) {
    
    curr.row = i + start.row - 1
    
    
    lapply(c(paste0("Starting analysis of sample ", MS.fastq.data$exp.name[curr.row], " at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # Go grab the relevant information about how the MIPS did 
    MIP.data[[i]] = get.MIP.Barcodes(MIP.info, MS.fastq.loc[curr.row], MS.fastq.val[curr.row])    
    
    # Now cbind it to the mapping information to save later
    mapping.results = cbind(mapping.results, unlist(MIP.data[[i]]$snv.wt.unique.ratio))
    
    # With each iteration we should try to save data about each MIP for each experiment
    # For now the number of reads would be good to know
    read.summary = cbind(read.summary, unlist(MIP.data[[i]]$total.match.lig))
    
    # Now save this table to the running directory      
    write.table(MIP.data[[i]], file=paste(c(output.prefix.list[curr.row], "results", "txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  }
  
  colnames(mapping.results) = c("MIP.name", "snv.strain", MS.fastq.data$exp.name[start.row:end.row])
  
  # Now collate all the data into a larger file of information
  all.MIP.info = summarize.MIPS(MIP.info, MIP.data)
  write.table(all.MIP.info, file=paste(c(output.directory, "MIP.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  write.table(mapping.results, file=paste(c(output.directory, "mappingdata.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  write.table(read.summary, file=paste(c(output.directory, "read.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  # At this point it would be nice to collate the data more into a single value summary per strain as well.
  # So take the average of the MIPs for a single strain and then generate a single value
  
  strain.average = get.strain.average(paste(c(output.directory, "mappingdata.txt"), collapse=""))
  write.table(strain.average[[1]], file=paste(c(output.directory, "strain.average.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  write.table(strain.average[[2]], file=paste(c(output.directory, "strain.average.sd.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)    

  lapply(c(paste0("Analysis complete at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  return (list(MIP.data, all.MIP.info, mapping.results))
}

get.strain.average = function(mapping.file) {
  
  # Read in the mapping file. It allows all the data to be properly converted to a numeric for the next steps
  mapping.results = read.table(mapping.file, header=TRUE, sep="\t")
  # Now replace all the "Na" with a 0 for the purposes of determining the average. 
  mapping.results[is.na(mapping.results)] = 0
  
  strain.list = unlist(unique(mapping.results$snv.strain))
  
  average.results = list()
  sd.results = list()
  
  for (i in 1:length(strain.list)) {
    
    curr.strain = strain.list[i]
    
    curr.rows = as.data.frame(mapping.results[mapping.results$snv.strain %in% curr.strain, 3:ncol(mapping.results)])
    
    average.results = rbind(average.results, apply(curr.rows, 2, mean))
    sd.results = rbind(sd.results, apply(curr.rows, 2, sd))
  }

  average.results = cbind(as.character(strain.list), average.results)
  sd.results = cbind(as.character(strain.list), sd.results)
  
  colnames(average.results) = c("snv.strain", colnames(mapping.results[3:ncol(mapping.results)]))
  colnames(sd.results) = c("snv.strain", colnames(mapping.results[3:ncol(mapping.results)]))
  
  average.results[is.na(average.results)] = 0
  sd.results[is.na(sd.results)] = 0
  
  return (list(average.results, sd.results))
  
}


### ----------------------- Functions related to the initial analysis and eventual formatting of data -----------------------###

# This function simply takes the given mapping positions (chromosome, locus) and returns a list of
# remapped positions along a single giant chromosome
generate.graph.data = function(name, chr, loci) {
  
  new.loci = list()
  
  for (i in 1:length(chr)) {
    
    if (chr[i] == "I") { new.loci[i] = loci[i] }
    else if (chr[i] == "II") { new.loci[i] = loci[i] + sum(LG.sizes[1]) }
    else if (chr[i] == "III") { new.loci[i] = loci[i] + sum(LG.sizes[1:2]) }
    else if (chr[i] == "IV") { new.loci[i] = loci[i] + sum(LG.sizes[1:3]) }
    else if (chr[i] == "V") { new.loci[i] = loci[i] + sum(LG.sizes[1:4]) }
    else if (chr[i] == "X") { new.loci[i] = loci[i] + sum(LG.sizes[1:5]) }
        
  }
  
  return (cbind(name, unlist(new.loci)))
}

generate.mapping.graph = function(mapping.data) {
  
  map.melt = melt(as.data.frame(mapping.data), id.vars="location")
  
  ggplot(map.melt, aes(x=location, y=value, colour=variable)) + geom_line()
  
  
  
}


# This takes in the MIP.data table and merely adds information to it for each MIP entry
# It begins by looking for the specific ligation arm sequence within a read
# Then it looks for gap fill reads matching the reference vs. mutant allele
# Then it checks to see how many of these have unique barcodes

# 140101 Update: To update the algorithm, we really need to make a change. It looks like with increasing MIPs,
# we may encounter issues where the same unique barcode has 1 or two reads, one of which is WT/SNV, the other having
# an additional mutation that causes it to be identified as a mismatch.
# 
# To remedy this let's try to identify unique barcodes BEFORE hand and compare just the sequence of the expected SNV?
# From examining some samples it appears that the barcodes remain unique amongst the mistmatches as well.

# 170222: Looking at data run, it appears that as we go over 1M reads per file, the length of time to run grows more than linearly
# 1M reads ~=20 minutes, 1.6M reads approaches 1H, 0.3M reads = 3minutes. Is there a way to speed up the process when working with huge data files?


get.MIP.Barcodes = function(MIP.data, MS.fastq.file, MS.fastq.values) {
  
  # Begin by taking the table and addin the appropriate new columns
  # We are most interested in how many total reads for any one allele, vs how many unique reads
  # We also want to gather allele ratios and determine how efficient this MIP is by how many reads overall match what we are looking for
  MIP.info = MIP.data
  # MIP.info$wt.bc.total = as.numeric(0)
  MIP.info$wt.bc.unique = as.numeric(0)
  # MIP.info$wt.PCR.ratio = as.numeric(0)
  
  #MIP.info$wt.bc.seq = ""
  
  # MIP.info$snv.bc.total = as.numeric(0)
  MIP.info$snv.bc.unique = as.numeric(0)
  # MIP.info$snv.PCR.ratio = as.numeric(0)
  
  #MIP.info$snv.bc.seq = ""
  
  # MIP.info$snv.wt.total.ratio = as.numeric(0)
  MIP.info$snv.wt.unique.ratio = as.numeric(0)
  
  
  MIP.info$total.match.lig = as.numeric(0)
  MIP.info$percent.mismatch = as.numeric(0)
  MIP.info$percent.pool = as.numeric(0)
  
  # 140818: Switch to identify high-quality matching percentage from original matching
  #MIP.info$low.qual.matching = as.numeric(0)
  #MIP.info$low.qual.percent = as.numeric(0)
  MIP.info$high.qual.match = as.numeric(0)
  MIP.info$high.qual.percent = as.numeric(0)
  # total.raw.bc is the number of unique high.qual.matche barcodes
  MIP.info$total.raw.bc = as.numeric(0)
  MIP.info$raw.PCR.ratio = as.numeric(0)
  
  
  MIP.info$non.dup.bc = as.numeric(0)  
  MIP.info$total.dup.bc = as.numeric(0)
  MIP.info$percent.dup.bc = as.numeric(0)
  MIP.info$unique.dup.bc = as.numeric(0)
  MIP.info$unique.dup.bc.accepted = as.numeric(0)
  MIP.info$total.unique.bc.accepted = as.numeric(0)
  MIP.info$percent.unique.bc.accepted = as.numeric(0)
  
  MIP.info$wt.bc.fuzzy = as.numeric(0)
  MIP.info$snv.bc.fuzzy = as.numeric(0)
  
  #wt.start.col = which(colnames(MIP.info) == "wt.bc.total")
  #wt.end.col = which(colnames(MIP.info) == "wt.PCR.ratio")
  #snv.start.col = which(colnames(MIP.info) == "snv.bc.total")
  #snv.end.col = which(colnames(MIP.info) == "snv.PCR.ratio")
  
  
  # Read in the given fastq sequence data file. This is held in memory for the entire analysis of the file!
  curr.fastq.seq = as.list(scan(file=MS.fastq.file, what=(seq="")))
  curr.fastq.val = as.list(scan(file=MS.fastq.values, what=(val="")))
  
  # Then you need to look at each MIP. For now just look at the total uniques for a single MIP allele
  # Layout of process:
  # 1. Pull down the matches by ligation sequence
  # 2. Filter poor quality reads at the SNV site
  # 3. Filter out duplicates and remove those in conflict
  # 4. Count matching WT (ref) reads, matching SNV reads
  
  for (i in 1:nrow(MIP.info)) {
    
    # Define your search items
    curr.lig.arm = toupper(MIP.info$lig.seq.read[i])
    curr.lig.length = nchar(curr.lig.arm)
    #lapply(c(paste0("curr.lig.arm: ", curr.lig.arm)), write, output.info, append=TRUE, ncolumns=1000)
    
    
    # If we are dealing with deletions, this information could invariably change. What we want is the correct length
    # Future proofed in case the read is shorter than expected due to deletions
    # Also good if MIP.read.length changes to say 75bp or more, then we only carry the gap.fill.seq that we care about
    curr.read.wt = toupper(MIP.info$wt.gap.fill.read[i])
    #lapply(c(paste0("curr.read.wt is: ", curr.read.wt)), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.read.wt = substr(curr.read.wt, 1, (MIP.read.length - curr.lig.length))
    #lapply(c(paste0("curr.read.wt after substr: ", curr.read.wt)), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.read.snv = toupper(MIP.info$snv.gap.fill.read[i])
    #lapply(c(paste0("curr.read.snv is: ", curr.read.snv)), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.read.snv = substr(curr.read.snv, 1, (MIP.read.length - curr.lig.length))
    #lapply(c(paste0("curr.read.snv after substr: ", curr.read.snv)), write, output.info, append=TRUE, ncolumns=1000)
    
    # WT lig + gap.fill read
    curr.MIP.wt = paste(c(curr.lig.arm, curr.read.wt), collapse="")
    # SNV lig + gap.fill read
    curr.MIP.snv = paste(c(curr.lig.arm, curr.read.snv), collapse="")
    
    # 1. Pull down the matches by ligation sequence
    
    # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
    # Create a map of matching sequences, then also grab their matching values
    # 150414: Changed the grep call to only look at the correct substring of bases for a match ie: 13-32 rather than searching all 50bp
    # Right now it may find the ligation arm at an incorrect position but then we'd have to check for matching subsets of the gap-fill sequence too.
    # I'd rather only count exact matches, and only use these as the correct number of reads to generate the "mismatch" percentage.
    # matching.lig.pos = grep(curr.lig.arm, curr.fastq.seq)
    
    # 160226: Should try to implement a fuzzy matching here where you can choose 1 or more bases to randomly mismatch except at the specific SNVs?
    #matching.list.pos = fuzzy.match(curr.lig.arm, do.call("substr", list(curr.fastq.seq, (MIP.barcode.length + 1), (MIP.barcode.length + curr.lig.length))), )
    
#     # We can probably make grep even faster if we slice the data down on each line to just position 13-42 (30bp that should contain the maximum ligation arm length)
#     curr.fastq.seq.data = unlist(lapply(curr.fastq.seq, function(x) substr(x, MIP.barcode.length+1, MIP.barcode.length + j)))
#     lapply(c(paste0("grep start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     
#     # See if we can grep for the data and get a list of positions back from it
#     # We don't really need any of the other original data for our purposes. We can just generate the seq and val files right now
#     # Use grepl to get a vector of TRUE/FALSE. We can then use this to kind of grab the other data... I think
#     # This will return all the positions of the matches with respect to the curr.fastq.seq.data set which has been trimmed in total entry length
#     # Note also that curr.fastq.seq is now truncated to only be the sequence data
#     
#     # Which positions match the exact grep call?
#     MIP.grep.match = grepl(curr.MIP.lig, curr.fastq.seq.data, perl=TRUE)
#     lapply(c(paste0("grep completed: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     
#     # Update the match list
#     MIP.match.loc = c(MIP.match.loc, which(MIP.grep.match == TRUE))
#     lapply(c(paste0("Memory usage for MIP.match.loc: ", object.size(MIP.match.loc))), write, output.info, append=TRUE, ncolumns=1000)
#     
    
        
    # Memory reduction call? the do.call returns a list that is super transient and should not continue to take up memory. But does it?
    matching.lig.pos = grep(curr.lig.arm, do.call("substr", list(curr.fastq.seq, (MIP.barcode.length + 1), (MIP.barcode.length + curr.lig.length))), perl=TRUE)
    matching.lig.seq = curr.fastq.seq[matching.lig.pos]
    matching.lig.val = curr.fastq.val[matching.lig.pos]

    # Garbage collect just in case that large list isn't purged from memory?
    gc()
    
    #     # Now we need to locate the specific basepair of interest. We'll evaluate it's quality before we even look at any other information
    #     curr.SNV.pos = MIP.barcode.length + curr.lig.length + as.numeric(MIP.info$lig.gap[i]) + 1
    #     # Grab only the SNV position that we're interested in terms of quality scores    
    #     matching.lig.val = unlist(lapply(matching.lig.val, function(x) substr(x, curr.SNV.pos, curr.SNV.pos)))
    #     
    #     # Now check them for quality
    #     matching.lig.pass.qual = lapply(matching.lig.val, function(x) all(get.quality.score(x)>=min.quality))
    #     # But now we won't know how many we're potentially losing from the total pool.
    #     matching.lig.seq = matching.lig.seq[unlist(matching.lig.pass.qual)] 
    

    # Subset into barcodes and reads - This could be memory intensive since each list is n elements long
    # Since the ligation arm could be 16-24 basepairs, we need to take it's length into account
    # The barcodes are always 12 basepairs long
    matching.lig.barcodes = unlist(lapply(matching.lig.seq, function (x) substr(x, 1, MIP.barcode.length)))
    #matching.lig.barcodes.val = unlist(lapply(matching.lig.val, function (x) substr(x, 1, MIP.barcode.length)))
            
    # Then we take from the end of the barcodes onwards
    matching.lig.reads = unlist(lapply(matching.lig.seq, function(x) substr(x, MIP.barcode.length + 1, MIP.read.length)))
    matching.lig.reads.val = unlist(lapply(matching.lig.val, function(x) substr(x, MIP.barcode.length + 1, MIP.read.length)))   
    
    # The total reads with matching ligation sequences
    MIP.info$total.match.lig[i] = length(matching.lig.seq)
    
    # 2. Filter poor quality reads at the SNV site
    # UPDATED 150914: Altered to evaluate multiple positions within the gap fill region. This is mainly to accomodate rDNA MIP analysis where our plasmid has multiple changes
    # Now lig.gap entry could be "x;y;z" format
    # UPDATED 160224: Looks like the location was incorrect for the values. We're identifying the gap value location and then looking at values
    # starting with 1 as the start of the ligation arm!! This won't work properly! We need to go in further to the correct location

    # ----- 160224: evaluate only gaps that have all bases > min.quality -----#
#     matching.lig.val = unlist(lapply(matching.lig.reads.val, function(x) substr(x, (MIP.barcode.length + curr.lig.length + 1), MIP.read.length)))
#     matching.lig.pass.qual = lapply(matching.lig.val, function(x) all(get.quality.score(x)>=min.quality))
#     # Now only examine these ones by updating the reads and barcodes
#     matching.lig.reads = matching.lig.reads[unlist(matching.lig.pass.qual)]
#     matching.lig.barcodes = matching.lig.barcodes[unlist(matching.lig.pass.qual)]
#     #matching.gap.seq = matching.gap.seq[!unlist(lapply(matching.gap.val, function(x) any(unlist(get.quality.score(x)) < 20)))]
    # ----- 160224: evaluate only gaps that have all bases > min.quality -----#
    
    # ----- Commented 160224 to see what evaluating all for q>=30 looks like -----    
    # Note we are adding the 1 here because normally if the entry is "0" then the very first base of the gap fill is the SNV site
    # 160224 changed this calculation to include the curr.lig.length!!! Super important if evaluating on quality of the actual SNV!!!
    
    lig.gap.val = as.numeric(unlist(strsplit(MIP.info$lig.gap[i], split=";"))) + 1 + curr.lig.length
    #lapply(c(paste0("lig.gap.val: ", lig.gap.val)), write, output.info, append=TRUE, ncolumns=1000)
    
    # Now we need to iterate through the possible locations and check them for value. This will definitely alter some values
    for (j in 1:length(lig.gap.val)) {
      
      # Now we need to locate the specific basepair of interest. We'll evaluate it's quality before we even look at any other information
      curr.SNV.pos = lig.gap.val[j]
      # Grab only the SNV position that we're interested in terms of quality scores    
      matching.lig.val = unlist(lapply(matching.lig.reads.val, function(x) substr(x, curr.SNV.pos, curr.SNV.pos)))
      # Now check them for quality
      matching.lig.pass.qual = lapply(matching.lig.val, function(x) all(get.quality.score(x)>=min.quality))
      
      # Now only examine these ones by updating the reads and barcodes
      matching.lig.reads = matching.lig.reads[unlist(matching.lig.pass.qual)]
      matching.lig.barcodes = matching.lig.barcodes[unlist(matching.lig.pass.qual)]
      
    }  
    
    # This large list is not used after this point. Destroy it. it should clean up the large matching.lig.pass.qual list as well.
    remove(matching.lig.reads.val)
    gc()
    
    # ----- end comment 160224 to see what evaluating all for q>=30 looks like -----    
    
    
    #     # Now we need to locate the specific basepair of interest. We'll evaluate it's quality before we even look at any other information
    #     curr.SNV.pos = curr.lig.length + as.numeric(MIP.info$lig.gap[i]) + 1
    #     # Grab only the SNV position that we're interested in terms of quality scores    
    #     matching.lig.val = unlist(lapply(matching.lig.reads.val, function(x) substr(x, curr.SNV.pos, curr.SNV.pos)))
    #     # Now check them for quality
    #     matching.lig.pass.qual = lapply(matching.lig.val, function(x) all(get.quality.score(x)>=min.quality))
    #     
    #     # Now only examine these ones by updating the reads and barcodes
    #     matching.lig.reads = matching.lig.reads[unlist(matching.lig.pass.qual)]
    #     matching.lig.barcodes = matching.lig.barcodes[unlist(matching.lig.pass.qual)]
    
    MIP.info$high.qual.match[i] = length(matching.lig.reads)
    MIP.info$high.qual.percent[i] = length(matching.lig.reads)/length(matching.lig.seq)
    MIP.info$total.raw.bc[i] = length(unlist(unique(matching.lig.barcodes)))
    MIP.info$raw.PCR.ratio[i] = as.numeric(MIP.info$high.qual.match[i])/as.numeric(MIP.info$total.raw.bc[i])
    
    # 3. Filter out duplicate barcodes and remove those in conflict
    
    # 140819 update ******************************************
    # Which barcodes are duplicated from the high-quality set
    # 160302 updated to 
    # 1) remove homopolymers of 4+ length within barcodes
    # 2) merge barcodes that are off by 1 mismatch?
    
    #1. Remove homopolymers
    matching.lig.barcodes.homopoly = unlist(lapply(matching.lig.barcodes, function(x) (check.homopolymers(x, max.poly) == 0)))
    # Update the barcode and reads
    matching.lig.barcodes = matching.lig.barcodes[matching.lig.barcodes.homopoly]
    matching.lig.reads = matching.lig.reads[matching.lig.barcodes.homopoly]
    
    #2. clean.MIP.barcodes should not only identify duplicates but also any 1-base mismatches
    # 160302: can't implement appropriately at this time. Just go with duplicate removal. At least we've removed homopolymers.
    #matching.lig.barcodes.dup.pos = clean.MIP.barcodes(matching.lig.barcodes)
    
    matching.lig.barcodes.dup.pos = duplicated(matching.lig.barcodes)

    # Grab the actual unique barcode sequences
    matching.lig.barcodes.dup.seq = unlist(unique(matching.lig.barcodes[matching.lig.barcodes.dup.pos]))
    #lapply(c(paste0("duplicated barcodes: ", head(matching.lig.barcodes.dup.seq))), write, output.info, append=TRUE, ncolumns=1000)
        
    # Now to update the list again to have ALL of the duplicated elements and not just the duplicated instances
    matching.lig.barcodes.dup.pos = (matching.lig.barcodes %in% matching.lig.barcodes.dup.seq)
    
    matching.lig.barcodes.dup.elements = matching.lig.barcodes[matching.lig.barcodes.dup.pos]
    matching.lig.reads.dup.elements = matching.lig.reads[matching.lig.barcodes.dup.pos]
    #lapply(c(paste0("duplicated barcodes elements: ", head(matching.lig.barcodes.dup.elements))), write, output.info, append=TRUE, ncolumns=1000)
    #lapply(c(paste0("duplicated barcodes lig reads: ", head(matching.lig.reads.dup.elements))), write, output.info, append=TRUE, ncolumns=1000)
    
        
    # Take all the unique barcode data as well
    matching.lig.barcodes.unique.elements = matching.lig.barcodes[!matching.lig.barcodes.dup.pos]
    matching.lig.reads.unique.elements = matching.lig.reads[!matching.lig.barcodes.dup.pos]
    
    MIP.info$non.dup.bc[i] = length(matching.lig.barcodes.unique.elements)
    
    # Now send this data to some kind of function to analyse separately and send back a list of two lists
    # 160224: This can't be working correctly if there is more than 1 SNV present for a MIP.
    # Fixed to resolve duplicates correctly
    matching.lig.unique.consensus = resolve.Duplicates(matching.lig.barcodes.dup.elements, matching.lig.reads.dup.elements, lig.gap.val)
    #matching.lig.unique.consensus = resolve.Duplicates(matching.lig.barcodes.dup.elements, matching.lig.reads.dup.elements, curr.SNV.pos)
    
    matching.lig.barcodes.unique.elements = c(matching.lig.barcodes.unique.elements, matching.lig.unique.consensus[[1]])
    matching.lig.reads.unique.elements = c(matching.lig.reads.unique.elements, matching.lig.unique.consensus[[2]])
    
    # Update the statistical information as we generate it 
    
    # total.dup.bc is the number of high qual duplicated barcode elements
    MIP.info$total.dup.bc[i] = length(matching.lig.barcodes.dup.elements)
    
    # 170410 adjusted this value to have the denominator as the number of raw.unique.barcodes
    MIP.info$percent.dup.bc[i] = as.numeric(MIP.info$total.dup.bc[i])/as.numeric(MIP.info$total.raw.bc[i])
    #MIP.info$percent.dup.bc[i] = as.numeric(MIP.info$total.dup.bc[i])/as.numeric(MIP.info$high.qual.match[i])
    MIP.info$unique.dup.bc[i] = length(matching.lig.barcodes.dup.seq)
    MIP.info$unique.dup.bc.accepted[i] = length(matching.lig.unique.consensus[[1]])
    MIP.info$total.unique.bc.accepted[i] = length(matching.lig.barcodes.unique.elements)
    MIP.info$percent.unique.bc.accepted[i] = as.numeric(MIP.info$total.unique.bc.accepted[i])/as.numeric(MIP.info$total.raw.bc[i])
    MIP.info$num.homopoly.bc[i] = as.numeric(length(matching.lig.barcodes.homopoly) - sum(matching.lig.barcodes.homopoly, na.rm=TRUE))
    
        
    #**********************************************************
    
    # 4. Count matching WT (ref) reads, matching SNV reads
    # Note that the final totals of wt.bc.unique and snv.bc.unique may not match those of total.unique.bc.accepted
    # At this point we've only filtered by matching the LIGATION arm sequence.
    # Quality was filtered by just looking at the position and not sequence of the SNV
    # Duplicates were self-compared so a lot of high-qual, non-matching gap sequences could slip through.
    # This is where $percent.mismatch helps us guage the quality of the MIP.
    
    # Now start filling the table with information      
    # MIP.info[i, wt.start.col:wt.end.col] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.MIP.wt)
    # MIP.info[i, snv.start.col:snv.end.col] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.MIP.snv)
    
    lig.gap.pos = as.numeric(unlist(strsplit(MIP.info$lig.gap[i], split=";"))) + 1
    
    # Generate the fuzzy.match versions of the wt and mut gap.fill seq
    #curr.read.wt.fuzzy = make.fuzzy.query(curr.read.wt, lig.gap.val, 1)
    #curr.read.snv.fuzzy = make.fuzzy.query(curr.read.snv, lig.gap.val, 1)

    # identifyMatches will already receive unique barcode information that's been duplicate.resolved and quality filtered at the SNv position(s).
    # It will not, however have looked at the specific sequences of the gap-fill reads.
    # Theoretically matching.lig.reads.unique.elements is composed of reads starting at the ligation arm through to the end of the read
    # We are sending lig.gap.pos based on it's location relative to the end of the ligation arm
    # 170410: FIXED this to use lig.gap.val instead which is relative position from the start of the read(start of lig.read)
    # curr.read.wt/snv is the gap-fill sequence
    
    # 170411 removed the need to specify the number of mismatches. Made a global variable instead that is set when the program is initiated.
    # This will allow us to change the fuzzy match off for regular mapping and population assessment.
    
    wt.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.wt, lig.gap.pos)
    snv.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.snv, lig.gap.pos)    
    
    #wt.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.wt, lig.gap.pos, 1)
    #snv.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.snv, lig.gap.pos, 1)
    
    #wt.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.wt, lig.gap.val, 1)
    #snv.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.snv, lig.gap.val, 1)
    
    MIP.info$wt.bc.unique[i] = wt.match.info[2]
    MIP.info$snv.bc.unique[i] = snv.match.info[2]
    
    MIP.info$wt.bc.fuzzy[i] = wt.match.info[3]
    MIP.info$snv.bc.fuzzy[i] = snv.match.info[3]
    
    #MIP.info$wt.bc.unique[i] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.wt, lig.gap.pos, 1)
    #MIP.info$snv.bc.unique[i] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.snv, lig.gap.pos, 1)
    
#     MIP.info$wt.bc.unique[i] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.MIP.wt, lig.gap.val, 1)
#     MIP.info$snv.bc.unique[i] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.MIP.snv, lig.gap.val, 1)
    
    
    # This was the total reads with matching ligation sequences
    # MIP.info$total.match.lig[i] = length(matching.lig.seq)
    # Updated to total reads with unique barcodes and matching ligation sequences that pass FIX THIS!
    # MIP.info$total.match.lig[i] = length(unlist(unique(matching.lig.barcodes)))
    
    #MIP.info$snv.wt.total.ratio[i] = as.numeric(MIP.info$snv.bc.total[i])/(as.numeric(MIP.info$snv.bc.total[i]) + as.numeric(MIP.info$wt.bc.total[i]))
    MIP.info$snv.wt.unique.ratio[i] = as.numeric(MIP.info$snv.bc.unique[i])/(as.numeric(MIP.info$snv.bc.unique[i]) + as.numeric(MIP.info$wt.bc.unique[i]))
  
    #MIP.info$dup.conflict.bc[i] = (length(matching.lig.barcodes.dup.seq) - length(matching.lig.unique.consensus[[1]]))
    #MIP.info$dup.conflict.bc.percent[i] = (MIP.info$dup.conflict.bc[i]/MIP.info$total.match.lig[i])
    
    #MIP.info$percent.mismatch[i] = (as.numeric(MIP.info$total.match.lig[i])- as.numeric(MIP.info$wt.bc.total[i]) - as.numeric(MIP.info$snv.bc.total[i]))/as.numeric(MIP.info$total.match.lig[i])
    MIP.info$percent.mismatch[i] = (as.numeric(MIP.info$total.unique.bc.accepted[i])- as.numeric(MIP.info$wt.bc.unique[i]) - as.numeric(MIP.info$snv.bc.unique[i]))/as.numeric(MIP.info$total.unique.bc.accepted[i])
    
    # This only gives us a number of reads from the total.match.lig that were removed by value. This would include mismatch instances as well
    # ************ We will also count read lost to SNV conflicts on duplicate barcodes as part of the low.qual.matching score *********
    #MIP.info$low.qual.matching[i] = MIP.info$total.match.lig[i]-length(matching.lig.reads)
    #MIP.info$low.qual.percent[i] = (MIP.info$low.qual.matching[i]/MIP.info$total.match.lig[i])
  
    
  }
  
  total.reads = sum(as.numeric(MIP.info$total.match.lig))
  
  # You can calculate what percentage of the experimental pool each MIP is taking up.
  MIP.info$percent.pool = sapply(MIP.info$total.match.lig, function(x) as.numeric(x)/total.reads)
  
  # Return the table of information
  return (MIP.info)
  
}

# 160302: New function to look find closely matching barcodes
# Could use the fuzzy mismatch algorithm to find partial duplicates BUT it become complicated
# AAGGTTCC -> AAGCTTCC -> AAGCTTCA
# The middle barcode is a match for both the first and third barcodes BUT 1 and 3 are now 2 mismatches from each other.
# How do I resolve this data? How do I choose which way the barcode goes? Do a complicated half-value scheme?
# Deprecate this idea until later.
clean.MIP.barcodes = function(barcode.data) {
  
  
  
}

# This is a new function (140819) used to resolve barcode duplicates and any potential conflicts
# 160224: Updated to handle read.SNV.pos with multiple positions
resolve.Duplicates = function(barcode.elements, read.elements, read.SNV.pos) {
  
  # Grab the list of repeated barcodes
  barcode.seq.list = unique(barcode.elements)
  new.barcode.elements = list()
  new.read.elements = list()

  # iterate through each barcode 
  for (i in 1:length(barcode.seq.list)) {
    
    matching.lig.SNV = list()
    curr.lig.reads = read.elements[barcode.elements %in% barcode.seq.list[i]]
    #lapply(c(paste0("head(curr.lig.reads): ", head(curr.lig.reads))), write, output.info, append=TRUE, ncolumns=1000)
    
        
    for (j in 1:length(read.SNV.pos)) {
      
      matching.lig.SNV = paste(matching.lig.SNV, unlist(lapply(curr.lig.reads, function(x) substr(x, read.SNV.pos[j], read.SNV.pos[j]))), sep="")
      
    }
    
    #lapply(c(paste0("head(matching.list.SNV): ", head(matching.lig.SNV))), write, output.info, append=TRUE, ncolumns=1000)
    
    # This generates a list of the possible SNVs read at this position
    # matching.lig.SNV = unlist(lapply(curr.lig.reads, function(x) substr(x, read.SNV.pos, read.SNV.pos)))
    # now turn it into a table so we can grab the mode
    SNV.table = table(as.vector(matching.lig.SNV))
    SNV.mode = names(SNV.table)[SNV.table == max(SNV.table)]
    
    # What if we have a conflict? ie same number of occurrences for each SNV? Drop the element completely otherwise pick the SNV.mode as the only read element for that barcode
    if (length(SNV.mode) == 1) {
      
      new.barcode.elements = c(new.barcode.elements, barcode.seq.list[i])
      
      # Take and read and replace the SNV with the mode SNV and use that as the new read element
      new.read = curr.lig.reads[1]
      substr(new.read, read.SNV.pos, read.SNV.pos) = SNV.mode
      new.read.elements = c(new.read.elements, new.read)      
    }    
    
  }
  
  # return the new duple to the caller
  return(list(new.barcode.elements, new.read.elements))
}

# This is primarily used to generate a summary file of the SNPs. It might not be all that useful in it's current state.
summarize.MIPS = function(MIP.info, MIP.data) {
  
  
  MIP.info.all = MIP.info
  #MIP.info.all$wt.bc.total = as.numeric(0)
  #MIP.info.all$wt.bc.sets = ""
  MIP.info.all$wt.bc.unique = as.numeric(0)
  MIP.info.all$wt.bc.unique.sets = ""
  #MIP.info.all$wt.bc.seq = ""
  
  #MIP.info.all$snv.bc.total = as.numeric(0)
  #MIP.info.all$snv.bc.sets = ""
  MIP.info.all$snv.bc.unique = as.numeric(0)
  MIP.info.all$snv.bc.unique.sets = ""
  #MIP.info.all$snv.bc.seq = ""    
  
  for (i in 1:length(MIP.data)) {
    
    for (j in 1:nrow(MIP.info.all)) {
      
      #MIP.info.all$wt.bc.total[j] = MIP.info.all$wt.bc.total[j] + as.numeric(MIP.data[[i]]$wt.bc.total[j])
      
      #MIP.info.all$wt.bc.sets[j] = paste(c(MIP.info.all$wt.bc.sets[j], MIP.data[[i]]$wt.bc.total[j]), collapse=";")
      
      MIP.info.all$wt.bc.unique[j] = MIP.info.all$wt.bc.unique[j] + as.numeric(MIP.data[[i]]$wt.bc.unique[j])
      
      MIP.info.all$wt.bc.unique.sets[j] = paste(c(MIP.info.all$wt.bc.unique.sets[j], MIP.data[[i]]$wt.bc.unique[j]), collapse=";")
      
      #MIP.info.all$wt.bc.seq[j] = paste(c(MIP.info.all$wt.bc.seq[j], MIP.data[[i]]$wt.bc.seq[j]), collapse="|")
      
      
      
      #MIP.info.all$snv.bc.total[j] = MIP.info.all$snv.bc.total[j] + as.numeric(MIP.data[[i]]$snv.bc.total[j])
      
      #MIP.info.all$snv.bc.sets[j] = paste(c(MIP.info.all$snv.bc.sets[j], MIP.data[[i]]$snv.bc.total[j]), collapse=";")
      
      MIP.info.all$snv.bc.unique[j] = MIP.info.all$snv.bc.unique[j] + as.numeric(MIP.data[[i]]$snv.bc.unique[j])
      
      MIP.info.all$snv.bc.unique.sets[j] = paste(c(MIP.info.all$snv.bc.unique.sets[j], MIP.data[[i]]$snv.bc.unique[j]), collapse=";")
      
      #MIP.info.all$snv.bc.seq[j] = paste(c(MIP.info.all$snv.bc.seq[j], MIP.data[[i]]$snv.bc.seq[j]), collapse="|")
      
    }
    
  } 
  
  return (MIP.info.all)
}


# Since we'll essentially treat WT and SNV sequences equally, we can run the same kind of algorithm to find information on them
# matching.lig.bardcodes: the list of 12bp barcodes that match up to the ligation reads
# matching.lig.reads: the list of 38bp reads that match up to the ligation reads
# curr.MIP.seq: the sequence we want to match against
# SNV.pos is based on the relative position AFTER the lig.arm (ie 1 = first base of curr.MIP.seq)

identifyMatches = function(matching.lig.barcodes, matching.lig.reads, curr.MIP.seq, SNV.pos){
  
  # 160226: Implement fuzzy matching on this grep function
  fuzzy.match.info = fuzzy.match(curr.MIP.seq, matching.lig.reads, SNV.pos, max.mismatch)
  matching.bc = matching.lig.barcodes[fuzzy.match.info[[1]]]
  
  #matching.bc = matching.lig.barcodes[fuzzy.match(curr.MIP.seq, matching.lig.reads, SNV.pos, num.mismatch)]
  
  #matching.bc = matching.lig.barcodes[grep(curr.MIP.seq, matching.lig.reads)]
  total.match = length(matching.bc)
  unique.match = length(unique(matching.bc))
  fuzzy.match = fuzzy.match.info[[3]]
  
  
  #results = c(length(matching.bc), length(unique(matching.bc)), paste(c(unique(matching.bc)), collapse=";"))
  #results = c(total.match, unique.match, total.match/unique.match )
  #results = total.match
  results = c(total.match, unique.match, fuzzy.match)
  
  return (results)
  
}

### ----------------------- Functions for normalizing the data after the initial analysis -----------------------###

# This function should take in a normalization file = which is basically a MIP.name x dilution table 
# $1 = MIP.name
# Need to somehow convert row names to a vector? Is that possible?
normalize.data = function(mapping.results.file, norm.file, output.file) {
  
  mapping.results = read.table(mapping.results.file, header=TRUE, sep="\t", colClasses = "character")
  new.mapping.results = mapping.results
  #norm.params = read.table(norm.file, header=FALSE, sep="\t", colClasses = "character")
  norm.params = read.table(norm.file, header=TRUE, sep="\t", colClasses = "character", row.names=1)
  type = "Mode"
  
  #   x.coord = norm.params[1, 2:ncol(norm.params)]
  #   norm.params = norm.params[2:nrow(norm.params), ]
  
  # Build the interpolation function for each MIP, then adjust all the values for that in the mapping.results information
  # 150121: Right now normalization isn't exclusively mapped to a specific MIP but rather by row. You should make it MIP-specific
  
  MIP.names = mapping.results$MIP.name
  
  for (i in 1:length(MIP.names)) {
    
    curr.MIP = MIP.names[i]
    mip.mean = as.numeric(norm.params[paste0(c(curr.MIP, "mean"), collapse="_"), type])
    mip.height = as.numeric(norm.params[paste0(c(curr.MIP, "h"), collapse="_"), type])
    mip.offset = as.numeric(norm.params[paste0(c(curr.MIP, "off"), collapse="_"), type])
    mip.sigma = as.numeric(norm.params[paste0(c(curr.MIP, "sigma"), collapse="_"), type])
    
    # Go through the table for that MIP and normalize all the values
    for (j in 3:ncol(mapping.results)) {
      
      new.mapping.results[i, j] = as.numeric(mapping.results[i,j]) - normalization.function(mip.mean, mip.sigma, mip.height, mip.offset, as.numeric(mapping.results[i, j]))
      #new.mapping.results[i, j] = normalization.function(mip.mean, mip.sigma, mip.height, mip.offset, as.numeric(mapping.results[i, j]))
    }    
    
    
  }
  
  
#   for (i in 1:nrow(mapping.results)) {
#     
#     #     mip.mean = as.numeric(norm.params$Mean[(2+4*(i-1))])
#     #     mip.height = as.numeric(norm.params$Mean[(3+4*(i-1))])
#     #     mip.offset = as.numeric(norm.params$Mean[(4+4*(i-1))])
#     #     mip.sigma = as.numeric(norm.params$Mean[(5+4*(i-1))])
#     
#     mip.mean = as.numeric(norm.params$Mode[(2+4*(i-1))])
#     mip.height = as.numeric(norm.params$Mode[(3+4*(i-1))])
#     mip.offset = as.numeric(norm.params$Mode[(4+4*(i-1))])
#     mip.sigma = as.numeric(norm.params$Mode[(5+4*(i-1))])
#     
#     # Go through the table for that MIP and normalize all the values
#     for (j in 3:ncol(mapping.results)) {
#       
#       new.mapping.results[i, j] = as.numeric(mapping.results[i,j]) - normalization.function(mip.mean, mip.sigma, mip.height, mip.offset, as.numeric(mapping.results[i, j]))
#       #new.mapping.results[i, j] = normalization.function(mip.mean, mip.sigma, mip.height, mip.offset, as.numeric(mapping.results[i, j]))
#     }    
#   }  
  
  write.table(new.mapping.results, file=output.file, append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  return (new.mapping.results)
  
}

normalization.function = function(mean, sigma, height, offset, x.val) {
  
  e = 2.71828183
  
  y.val = offset + height*e^(-0.5*(((x.val-mean)/sigma)^2))
  
  return (y.val)
  
}

### ----------------------- Sequence splitting files. Could potentially update or change -----------------------###

# Split files into a sequence vs quality file
# Techincally when Owen initially works with the fasta files, any failing clusters are not included BUT
# sequences with sporadic low values could still make it in through the process
# At this point we aren't investigating the quality value but it may be useful to look up later
# 141218: Need to fix this to be more efficient and line 2 = sequence, line 4 = values
# 150403: Found out that the scan function can read in .gz files! No more pre-unzipping them. Just use as is.
make.seq.files = function(MS.fastq.file, seq.directory, value.directory, start.row, end.row) {
  
  seq.directory = "seqFiles"
  value.directory = "seqValues"
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location
  #output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(directory, x), collapse="")))
  
  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], seq.directory), collapse="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], value.directory), collapse="/"))

  for (i in start.row:end.row) {
    
    #curr.fastq = scan(file=MS.fastq.loc[i], what=(seq=""), sep="\n")
    #curr.fastq.seq = curr.fastq[seq(2L, length(curr.fastq), 4L)] # grab lines starting at 2 and every 4th after that
    #curr.fastq.values = curr.fastq[seq(4L, length(curr.fastq), 4L)]

    split.size = split.seq.file(MS.fastq.loc[i], MS.fastq.data$seqLocation[i], MS.fastq.data$valLocation[i])
    
    MS.fastq.data$seq.file.size[i] = split.size
    MS.fastq.data$val.file.size[i] = split.size
    
    lapply(c(paste0(MS.fastq.loc[i], " reads:\t", split.size)), write, output.info, append=TRUE, ncolumns=1000)
    
    #write.table(curr.fastq.seq, file=MS.fastq.data$seqLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)      
    #write.table(curr.fastq.values, file=MS.fastq.data$valLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
  }

  write.table(MS.fastq.data, file=MS.fastq.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

}


# 170209 new version of make.seq.files based on faster opening of .gz files and writing of them. Should also reduce memory usage by a lot
split.seq.file = function(input.file, seq.file, val.file) {
  
  #lapply(c(paste0("read.fastq.gz test iteration start at: ", Sys.time())), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  #lapply(c(paste0("file to read: ", input.file)), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  #read.length = c(10000, 100000, 1000000)
  #read.length = c(1000000)
  
  #for (i in 1:length(read.length)) {
  
  nlines = 1000000
  
  input.con = file(input.file, "r")
  seq.con = gzfile(seq.file, "w")
  val.con = gzfile(val.file, "w")
  
  seq.return = character()
  val.return = character()
  file.size = 0
  
  #lapply(c(paste0("nlines: ", nlines)), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  #lapply(c(paste0(Sys.time(), ": start")), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  
  #while (length(curr.data <- scan(file=input.con, what=character(), nlines=4, sep="\n", quiet=TRUE)) > 0) {
  
  # Use >3 in this argument because a well-formed fastq file should always have 0 or at least 4 lines of data at the end. Anything else would be suspicious
  while (length(curr.data <- readLines(con=input.con, n=nlines)) > 3) {
    
    file.size = file.size + length(curr.data)/4
    
    seq.length = min(nlines, length(curr.data))
    
    #seq.return = c(seq.return, curr.data[seq(2, seq.length, 4)])
    #val.return = c(val.return, curr.data[seq(4, seq.length, 4)])
    
    cat(curr.data[seq(2, seq.length, 4)], file=seq.con, sep="\n")
    cat(curr.data[seq(4, seq.length, 4)], file=val.con, sep="\n")

  
  }
  
  #lapply(c(paste0(Sys.time(), ": end")), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  #lapply(c(paste0("Lines in file: ", length(seq.return))), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  #write.table(seq.return, file=paste0(input.file, ".read.fastq.", nlines, ".lines"), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
  close (input.con)
  close (seq.con)
  close (val.con)
  #return (list(seq.return, val.return))
  
  #}
  
  return (file.size)
}

# 161031: Used exclusively when working with Data mixed with RNAseq runs from Adam's NextSeq data.
# The main issue is that files are extremely large and we should try to pull out only MIP-matching reads before
# further generating the seq and val version of the files. 
# Steps:
# 1. Open up the files
# 2. Index lines where MIP ligation arms are matching
# 3. generate 4 files: MIP.data and RNA.data vs seq.data, val.data
# 4. save the files.
# make.MIP.seq.files("./MiSeq.analysis/Analysis.Run.files/161025.M20.NS.RNAseq.fastq.files.txt", "mix16.MIP.info.txt", "seqFiles", "seqValues", 1, 1)

make.MIP.seq.files = function(MS.fastq.file, mapping.info.file, seq.directory, value.directory, start.row, end.row) {
    
  lapply(c(paste0("Starting RNAseq demux: ", Sys.time())), write, "161104.RNAseq.demux.test.info", append=TRUE, ncolumns=1000)

  
  seq.directory = "seqFiles"
  value.directory = "seqValues"
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location
  #output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(directory, x), collapse="")))
  
  # Read the mapping information into a table examination
  MIP.info = read.table(mapping.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], seq.directory), collapse="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], value.directory), collapse="/"))
  
  # Make the search pattern from all the ligation arms (in upper case)
  # 161102: fixed to use the lig.seq.read column! It was using lig.seq previously and returning no real results.
  MIP.lig = paste(unlist(toupper(MIP.info$lig.seq.read)), collapse="|")
  
  #grep(paste(unlist(toupper(MIP.data$lig.seq)), collapse="|"), c("CTACTTACCCCATATGTTGTCG", "CTACTTACCCCATATGTTGTCG"))
  
  for (i in start.row:end.row) {
    
    # Scan the file in
    lapply(c(paste0("Scanning in current file: ", MS.fastq.loc[i])), write, "161104.RNAseq.demux.test.info", append=TRUE, ncolumns=1000)
    lapply(c(paste0("Scan start: ", Sys.time())), write, "161104.RNAseq.demux.test.info", append=TRUE, ncolumns=1000)
    curr.fastq = scan(file=MS.fastq.loc[i], what=(seq=""), sep="\n")
    
    lapply(c(paste0("Scan finish: ", Sys.time())), write, "161104.RNAseq.demux.test.info", append=TRUE, ncolumns=1000)
    
    # See if we can grep for the data and get a list of positions back from it
    # We don't really need any of the other original data for our purposes. We can just generate the seq and val files right now
    # Use grepl to get a vector of TRUE/FALSE. We can then use this to kind of grab the other data... I think
    MIP.match.loc = grepl(MIP.lig, curr.fastq)
    lapply(c(paste0("grep completed: ", Sys.time())), write, "161104.RNAseq.demux.test.info", append=TRUE, ncolumns=1000)
    
    MIP.fastq.seq.loc = which(MIP.match.loc == TRUE)
    MIP.fastq.val.loc = MIP.fastq.seq.loc + 2
    
    # Also rebuild the fastq file that contains the data you didn't want
    # grab the positions of the false sequence hits
    other.fastq.seq.loc = which(MIP.match.loc[seq(2, length(MIP.match.loc), 4)] == FALSE)
    other.fastq.info.loc = other.fastq.seq.loc - 1
    #other.fastq.plus.loc = other.fastq.seq.loc + 1
    #other.fastq.val.loc = other.fastq.seq.loc + 2
    
    other.fastq.all.loc = unlist(lapply(other.fastq.info.loc, function(x) seq(x, x+3, 1)))
    
    curr.fastq.seq = curr.fastq[MIP.fastq.seq.loc] # grab lines starting at 2 and every 4th after that
    curr.fastq.values = curr.fastq[MIP.fastq.val.loc]
    
    MS.fastq.data$seq.file.size[i] = length(curr.fastq.seq)
    MS.fastq.data$val.file.size[i] = length(curr.fastq.values)
    
    write.table(curr.fastq.seq, file=MS.fastq.data$seqLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)      
    write.table(curr.fastq.values, file=MS.fastq.data$valLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(curr.fastq[other.fastq.all.loc], file=paste(c(MS.fastq.loc[i], "demux.data.fastq"), sep=""), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
  }
  
  write.table(MS.fastq.data, file=MS.fastq.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}


# Currently the scan on a large file is maybe 5 minutes or so but the grep and file completion takes HOURS! How can we speed it up?
# 1. We don't need to read each line. Just read the sequencing data lines should speed it up 4X

make.MIP.seq.files.fast = function(MS.fastq.file, mapping.info.file, seq.directory, value.directory, start.row, end.row) {
  
  #library(data.table)
  # consider fread but can't open zipped files so abandoned
  output.info = "161104.RNAseq.demux.test.info"
  
  lapply(c(paste0("\nStarting RNAseq demux rerun: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  seq.directory = "seqFiles"
  value.directory = "seqValues"
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.R1loc = MS.fastq.data$R1location
  #MS.fastq.R2loc = MS.fastq.data$R2location
  #output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(directory, x), collapse="")))
  
  # Read the mapping information into a table examination
  MIP.info = read.table(mapping.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  dir.info = unlist(strsplit(MS.fastq.R1loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], seq.directory), collapse="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], value.directory), collapse="/"))
  
  # Make the search pattern from all the ligation arms (in upper case)
  # 161102: fixed to use the lig.seq.read column! It was using lig.seq previously and returning no real results.
  #MIP.lig = paste(unlist(toupper(MIP.info$lig.seq.read)), collapse="|")
  #MIP.ext = paste(unlist(toupper(MIP.info$ext.seq.read)), collapse="|")
  
  MIP.lig = toupper(MIP.info$lig.seq.read)
  MIP.ext = toupper(MIP.info$ext.seq.read)
  
  #grep(paste(unlist(toupper(MIP.data$lig.seq)), collapse="|"), c("CTACTTACCCCATATGTTGTCG", "CTACTTACCCCATATGTTGTCG"))
  
  for (i in start.row:end.row) {
    
    # Scan the file in
    lapply(c(paste0("Scanning in current file: ", MS.fastq.R1loc[i])), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Scan start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    curr.fastq = scan(file=MS.fastq.R1loc[i], what=(seq=""), sep="\n")
    
    lapply(c(paste0("Scan finish: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # Grab only the sequence data lines for grep
    curr.fastq.seq.loc = seq(2, length(curr.fastq), 4)
    lapply(c(paste0("grep list size: ", length(curr.fastq.seq.loc))), write, output.info, append=TRUE, ncolumns=1000)
    # Implement a memory saving measure by dropping all the unneeded data
    # curr.fastq.seq.data = curr.fastq[curr.fastq.seq.loc]
    curr.fastq.seq = curr.fastq[curr.fastq.seq.loc]
    curr.fastq.val = curr.fastq[(curr.fastq.seq.loc+2)]

    lapply(c(paste0("Memory usage for curr.fastq: ", object.size(curr.fastq))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for curr.fastq.seq: ", object.size(curr.fastq.seq))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for curr.fastq.val: ", object.size(curr.fastq.val))), write, output.info, append=TRUE, ncolumns=1000)
    
    remove(curr.fastq)
    
    # 161109 Update to do a fixed=TRUE kind of grep but this works only strictly (and much faster supposedly) on exact matches.
    # We must therefore generate a for loop (between 19 and 22 bp) of the sequences in order to capture all the possibilities
    # Can we reduce the number of samples to look through at each iteration? May require more work than it's worth since we are usually just a small
    # percentage of the data
    
    MIP.match.loc = integer()
    
    for (j in min(nchar(MIP.lig)):max(nchar(MIP.lig))) {
      
      lapply(c(paste0("grep attempt on ligation length: ", j)), write, output.info, append=TRUE, ncolumns=1000)
      
      # Which ligation sequences are of length J?
      #MIP.lig = paste(unlist(toupper(MIP.info$lig.seq.read)), collapse="|")
      
      curr.MIP.lig = paste(unlist(MIP.lig[which(nchar(MIP.lig) == j)]), collapse="|")
      
      # We can probably make grep even faster if we slice the data down on each line to just position 13-42 (30bp that should contain the maximum ligation arm length)
      curr.fastq.seq.data = unlist(lapply(curr.fastq.seq, function(x) substr(x, MIP.barcode.length+1, MIP.barcode.length + j)))
      lapply(c(paste0("grep start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
      
      # See if we can grep for the data and get a list of positions back from it
      # We don't really need any of the other original data for our purposes. We can just generate the seq and val files right now
      # Use grepl to get a vector of TRUE/FALSE. We can then use this to kind of grab the other data... I think
      # This will return all the positions of the matches with respect to the curr.fastq.seq.data set which has been trimmed in total entry length
      # Note also that curr.fastq.seq is now truncated to only be the sequence data
      
      # Which positions match the exact grep call?
      MIP.grep.match = grepl(curr.MIP.lig, curr.fastq.seq.data, perl=TRUE)
      lapply(c(paste0("grep completed: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
      
      # Update the match list
      MIP.match.loc = c(MIP.match.loc, which(MIP.grep.match == TRUE))
      lapply(c(paste0("Memory usage for MIP.match.loc: ", object.size(MIP.match.loc))), write, output.info, append=TRUE, ncolumns=1000)
      
      
      #MIP.fastq.seq.loc = curr.fastq.seq.loc[MIP.match.loc == TRUE]
      #MIP.fastq.val.loc = MIP.fastq.seq.loc + 2
      
    }
    
    # The matches are now all updated in the list and can be identified and in turn, we can identify those NOT in the matching list
    lapply(c(paste0("Match locations calculated: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    other.fastq.info.R1loc = curr.fastq.seq.loc[!curr.fastq.seq.loc %in% MIP.match.loc] - 1 
    
    
    # Also rebuild the fastq file that contains the data you didn't want
    # grab the positions of the false sequence hits
    #other.fastq.seq.loc = curr.fastq.seq.loc[MIP.match.loc == FALSE]
    #other.fastq.info.R1loc = other.fastq.seq.loc - 1
    #remove (other.fastq.seq.loc)
    lapply(c(paste0("Non-match fastq start locations indexed : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for other.fastq.info.R1loc: ", object.size(other.fastq.info.R1loc))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Matching lines: ", length(MIP.match.loc))), write, output.info, append=TRUE, ncolumns=1000)
    
    
    
    #other.fastq.all.loc = unlist(lapply(other.fastq.info.loc, function(x) seq(x, x+3, 1)))
    #lapply(c(paste0("Non-match all lines calculated : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    #curr.fastq.seq = curr.fastq[MIP.fastq.seq.loc] # grab lines starting at 2 and every 4th after that
    #curr.fastq.values = curr.fastq[MIP.fastq.val.loc]
    
    #MS.fastq.data$seq.file.size[i] = length(curr.fastq.seq)
    MS.fastq.data$seq.file.size[i] = length(MIP.match.loc)
    MS.fastq.data$val.file.size[i] = MS.fastq.data$seq.file.size[i]
    
    write.table(curr.fastq.seq[MIP.match.loc == TRUE], file=MS.fastq.data$seqLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    lapply(c(paste0("seq.fastq written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    write.table(curr.fastq.val[MIP.match.loc == TRUE], file=MS.fastq.data$valLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    lapply(c(paste0("val.fastq written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)

    
#     # Step 2: We need to redo the entire thing on the read2 data and compare/merge. This will be rather costly time-wise
#     # 1. remove the current fastq to clear up space and the other fastq files you just made because adding to them from the read2 data wouldn't give any more information
#     remove(curr.fastq, curr.fastq.seq, curr.fastq.values)
#     
#     # Scan the file in
#     lapply(c(paste0("Scanning in current R2 file: ", MS.fastq.R2loc[i])), write, output.info, append=TRUE, ncolumns=1000)
#     lapply(c(paste0("Scan start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     curr.fastq = scan(file=MS.fastq.R2loc[i], what=(seq=""), sep="\n")
#     
#     lapply(c(paste0("Scan finish: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     
#     # Grab only the sequence data lines for grep
#     curr.fastq.seq.loc = seq(2, length(curr.fastq), 4)
#     curr.fastq.seq.data = curr.fastq[curr.fastq.seq.loc]
#     
#     # We can probably make grep even faster if we slice the data down on each line to just position 13-42 (30bp that should contain the maximum ligation arm length)
#     #matching.lig.pos = grep(curr.lig.arm, do.call("substr", list(curr.fastq.seq, (MIP.barcode.length + 1), (MIP.barcode.length + curr.lig.length))))
#     curr.fastq.seq.data = unlist(lapply(curr.fastq.seq.data, function(x) substr(x, MIP.barcode.length+1, 42)))
#     
#     lapply(c(paste0("grep list size: ", length(curr.fastq.seq.data))), write, output.info, append=TRUE, ncolumns=1000)
#     lapply(c(paste0("grep start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     
#     # See if we can grep for the data and get a list of positions back from it
#     # We don't really need any of the other original data for our purposes. We can just generate the seq and val files right now
#     # Use grepl to get a vector of TRUE/FALSE. We can then use this to kind of grab the other data... I think
#     
#     MIP.match.loc = grepl(MIP.lig, curr.fastq.seq.data)
#     lapply(c(paste0("grep completed: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     
#     # Also rebuild the fastq file that contains the data you didn't want
#     # grab the positions of the false sequence hits
#     #other.fastq.seq.loc = curr.fastq.seq.loc[MIP.match.loc == FALSE]
#     other.fastq.info.R2loc = curr.fastq.seq.loc[MIP.match.loc == FALSE] - 1
#     lapply(c(paste0("Non-match start locations calculated : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     
#     lapply(c(paste0("Keep only non-match locations from both sets of data : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
#     lapply(c(paste0("R1 non-match : ", length(other.fastq.info.R1loc))), write, output.info, append=TRUE, ncolumns=1000)
#     lapply(c(paste0("R2 non-match : ", length(other.fastq.info.R2loc))), write, output.info, append=TRUE, ncolumns=1000)
#     other.fastq.info.all.loc = other.fastq.info.R1loc %in% other.fastq.info.R2loc
#     
#     write.table(other.fastq.info.all.loc, file=paste0(c(MS.fastq.R1loc[i], "dual.demux.data.loc"), collapse="."), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)    
#     #write.table(curr.fastq[other.fastq.all.loc], file=paste0(c(MS.fastq.R1loc[i], "demux.data.fastq"), collapse="."), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
#     lapply(c(paste0("RNAseq.file written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
  }
  
  write.table(MS.fastq.data, file=MS.fastq.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}

make.MIP.seq.files.fast.perl = function(MS.fastq.file, mapping.info.file, seq.directory, value.directory, start.row, end.row) {
  
  #library(data.table)
  # consider fread but can't open zipped files so abandoned
  output.info = "161110.RNAseq.demux.test.info"
  
  lapply(c(paste0("\nStarting RNAseq demux perl rerun: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  seq.directory = "seqFiles"
  value.directory = "seqValues"

  # We need to grab some data about file locations  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.R1loc = MS.fastq.data$R1location
  MS.fastq.R2loc = MS.fastq.data$R2location
  MS.fastq.R1demuxloc = MS.fastq.data$R1demux
  MS.fastq.R2demuxloc = MS.fastq.data$R2demux

  # Read the mapping information into a table examination
  MIP.info = read.table(mapping.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  dir.info = unlist(strsplit(MS.fastq.R1loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], seq.directory), collapse="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], value.directory), collapse="/"))
  
  # Make the search pattern from all the ligation arms (in upper case)
  # 161102: fixed to use the lig.seq.read column! It was using lig.seq previously and returning no real results.
  MIP.lig = paste(unlist(toupper(MIP.info$lig.seq.read)), collapse="|")
  MIP.ext = paste(unlist(toupper(MIP.info$ext.seq)), collapse="|")
  

  for (i in start.row:end.row) {
    
    # Scan the file in. This is the longest step and you should consider trying to make it better. 
    # It's also the most memory intensive step.
    lapply(c(paste0("Scanning in current file: ", MS.fastq.R1loc[i])), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Scan start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.fastq = scan(file=MS.fastq.R1loc[i], what=(seq=""), sep="\n")
    
    lapply(c(paste0("Scan finish: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # Grab only the sequence data lines for grep
    curr.fastq.seq.loc = seq(2, length(curr.fastq), 4)
    lapply(c(paste0("grep list size: ", length(curr.fastq.seq.loc))), write, output.info, append=TRUE, ncolumns=1000)
    
    # Implement a memory saving measure by dropping all the unneeded data
    
    curr.fastq.seq = curr.fastq[curr.fastq.seq.loc]
    curr.fastq.val = curr.fastq[(curr.fastq.seq.loc+2)]
    
    lapply(c(paste0("Memory usage for curr.fastq: ", object.size(curr.fastq))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for curr.fastq.seq: ", object.size(curr.fastq.seq))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for curr.fastq.val: ", object.size(curr.fastq.val))), write, output.info, append=TRUE, ncolumns=1000)
    
    remove(curr.fastq)
    gc()
    
    # We can probably make grep even faster if we slice the data down on each line to just position 13-42 (30bp that should contain the maximum ligation arm length)
    # This step is best accomplished with a straight call to substr since it's a list. It's the fastest by a gigantic factor based on testing.
    curr.fastq.seq.data = substr(curr.fastq.seq, MIP.barcode.length + 1, MIP.barcode.length + max(nchar(MIP.info$lig.seq.read)))
    
    # See if we can grep for the data and get a list of positions back from it
    # We don't really need any of the other original data for our purposes. We can just generate the seq and val files right now
    # Use grepl to get a vector of TRUE/FALSE. We can then use this to kind of grab the other data... I think
    # This will return all the positions of the matches with respect to the curr.fastq.seq.data set which has been trimmed in total entry length
    # Note also that curr.fastq.seq is now truncated to only be the sequence data
    
    lapply(c(paste0("grep start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # Which positions match the exact grep call?
    MIP.match.loc = grepl(MIP.lig, curr.fastq.seq.data, perl=TRUE)
    
    lapply(c(paste0("grep completed: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for MIP.match.loc: ", object.size(MIP.match.loc))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Matching lines: ", length(which(MIP.match.loc == TRUE)))), write, output.info, append=TRUE, ncolumns=1000)
        
    # The matches are now all updated in the list and can be identified and in turn, we can identify those NOT in the matching list
    lapply(c(paste0("Match locations calculated: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)

    # Update the fastq.size data
    MS.fastq.data$seq.file.size[i] = length(MIP.match.loc)
    MS.fastq.data$val.file.size[i] = MS.fastq.data$seq.file.size[i]
    
    # At this point you can write the seq and val files to disk. Looking at the 2nd read data won't change anything since the ligation mapping is so strict.
    write.table(curr.fastq.seq[MIP.match.loc], file=MS.fastq.data$seqLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    lapply(c(paste0("seq.fastq written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    write.table(curr.fastq.val[MIP.match.loc], file=MS.fastq.data$valLocation[i], append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    lapply(c(paste0("val.fastq written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # clean up the memory a bit and garbage collect
    remove(curr.fastq.seq, curr.fastq.val)
    gc()

    # ---------- Now to deal with pulling out the non-MIP data ---------- #
    
    # Index all the start locations of any non-matching data
    other.fastq.info.R1loc = curr.fastq.seq.loc[!MIP.match.loc] - 1
    remove (MIP.match.loc)

    lapply(c(paste0(length(other.fastq.info.R1loc), " R1 non-match fastq start locations indexed : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for other.fastq.info.R1loc: ", object.size(other.fastq.info.R1loc))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Proceed to 2nd read analysis at : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.fastq = scan(file=MS.fastq.R2loc[i], what=(seq=""), sep="\n")
    
    #curr.fastq.seq = curr.fastq[curr.fastq.seq.loc]
    # note that curr.fastq.seq.loc is already generated and was never removed so we can re-use that list
    curr.fastq.seq.data = substr(curr.fastq[curr.fastq.seq.loc], 1, max(nchar(MIP.info$ext.seq.read)))
    
    lapply(c(paste0("grep start: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # Which positions match the exact grep call?
    MIP.match.loc = grepl(MIP.ext, curr.fastq.seq.data, perl=TRUE)
    
    lapply(c(paste0("grep completed: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Memory usage for MIP.match.loc: ", object.size(MIP.match.loc))), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0("Matching lines: ", length(which(MIP.match.loc == TRUE)))), write, output.info, append=TRUE, ncolumns=1000)
    
    # The matches are now all updated in the list and can be identified and in turn, we can identify those NOT in the matching list
    other.fastq.info.R2loc = curr.fastq.seq.loc[!MIP.match.loc] - 1
    lapply(c(paste0(length(other.fastq.info.R2loc), " R2 non-match fastq start locations indexed : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    # Now compare with the R1 locations and only keep those that match
    other.fastq.final.loc = other.fastq.info.R1loc[other.fastq.info.R1loc %in% other.fastq.info.R2loc]
    lapply(c(paste0(length(other.fastq.final.loc), " Merged non-match fastq start locations indexed : ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    write.table(other.fastq.final.loc, file=paste0(c(MS.fastq.R1demuxloc[i], "dual.demux.data.loc"), collapse="."), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)    
    #write.table(curr.fastq[other.fastq.all.loc], file=paste0(c(MS.fastq.R1loc[i], "demux.data.fastq"), collapse="."), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    lapply(c(paste0("RNAseq.file written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)

    # ---------- Regenerate final fastq file ---------- #
    # At this point we should be able to regenerate the final fastq files, hopefully as .gx files using gzfile as a command?
    # Still have R2 fastq open so let's rewrite that first, and then re-open and re-write the R1 file
    # We have all the index start locations, now we need to expand that
    
    remove(curr.fastq)
    gc()
    
    write.fastq.gz(MS.fastq.R1loc[i], paste0(c(MS.fastq.R1demuxloc[i], "demux.data.fastq.gz"), collapse="."), other.fastq.final.loc)
    write.fastq.gz(MS.fastq.R2loc[i], paste0(c(MS.fastq.R2demuxloc[i], "demux.data.fastq.gz"), collapse="."), other.fastq.final.loc)
    
      
    #------------- altered this with the new pipeline ------------------#
    
    other.fastq.final.loc = unlist(lapply(other.fastq.final.loc, function(x) seq(x, x+3, 1)))
    
    # Open a connection to a new .gz file
    fastq.gz = gzfile(paste0(c(MS.fastq.R2demuxloc[i], "demux.data.fastq.gz"), collapse="."), "w")
    cat(curr.fastq[other.fastq.final.loc], file=fastq.gz, sep="\n")
    close(fastq.gz)
    lapply(c(paste0("R2 demuxed fastq file written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    remove(curr.fastq)
    gc()
    
    # Repeat for R1
    curr.fastq = scan(file=MS.fastq.R1loc[i], what=(seq=""), sep="\n")
    fastq.gz = gzfile(paste0(c(MS.fastq.R1demuxloc[i], "demux.data.fastq.gz"), collapse="."), "w")
    cat(curr.fastq[other.fastq.final.loc], file=fastq.gz, sep="\n")
    close(fastq.gz)
    lapply(c(paste0("R1 demuxed fastq file written: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
    
    remove(curr.fastq, other.fastq.final.loc)    
    
    #------------- altered this with the new pipeline ------------------#
    
    
    gc()
  }
  
  write.table(MS.fastq.data, file=MS.fastq.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}

# Designed specifically for faster fastq scan?, instead of reading in the entire file and then grabbing what we want, just do it from the start?
# Really just saves on initial memory needs for opening the file
# returns a seq.data and val.data set of lists
# relevant data is at line 2+4n (seq) and 4+4n (val)
# 161121: testing suggest we can read 10^6 lines at a time without any big issues but test file was less than 10^6 lines
#         Tried is on a 46M line fastq.gz file and it ran in just under 3 minutes!
read.fastq.gz = function(input.file) {
  
  lapply(c(paste0("read.fastq.gz test iteration start at: ", Sys.time())), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  lapply(c(paste0("file to read: ", input.file)), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  #read.length = c(10000, 100000, 1000000)
  read.length = c(1000000)
  
  #for (i in 1:length(read.length)) {
  
  nlines = read.length[i]
  input.con = file(input.file, "r")
  seq.return = character()
  val.return = character()
  
  
  #lapply(c(paste0("nlines: ", nlines)), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  #lapply(c(paste0(Sys.time(), ": start")), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  
  #while (length(curr.data <- scan(file=input.con, what=character(), nlines=4, sep="\n", quiet=TRUE)) > 0) {
  
  # Use >3 in this argument because a well-formed fastq file should always have 0 or at least 4 lines of data at the end. Anything else would be suspicious
  while (length(curr.data <- readLines(con=input.con, n=nlines)) > 3) {
    
    seq.length = min(nlines, length(curr.data))
    
    seq.return = c(seq.return, curr.data[seq(2, seq.length, 4)])
    val.return = c(val.return, curr.data[seq(4, seq.length, 4)])
    
  }
  
  #lapply(c(paste0(Sys.time(), ": end")), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  #lapply(c(paste0("Lines in file: ", length(seq.return))), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  #write.table(seq.return, file=paste0(input.file, ".read.fastq.", nlines, ".lines"), append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
  close (input.con)
  #return (list(seq.return, val.return))
  
  #}
  
  return (list(seq.return, val.return))
}

write.new.gz = function(data.info.file, output.info) {
  
  lapply(c(paste0("\nWrite.new.gz started at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  data.info.df = read.table(data.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  #demux.file.loc = data.info.df$demux.loc.file
  #R1.input = data.info.df$R1location
  #R2.input = data.info.df$R2location
  
  for (i in 1:nrow(data.info.df)) {
      
    demux.loc = unlist(scan(file=data.info.df$demux.loc.file[i], what=numeric()))

    write.fastq.gz(data.info.df$R1location[i], paste0(c(data.info.df$demux.loc.file[i], "161129.R1.demux.data.fastq.gz"), collapse="."), demux.loc, output.info)    
    write.fastq.gz(data.info.df$R2location[i], paste0(c(data.info.df$demux.loc.file[i], "161129.R2.demux.data.fastq.gz"), collapse="."), demux.loc, output.info)
    
    lapply(c(paste0("Re-write for set completed at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  }
  
  lapply(c(paste0("All data completed at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  
}


# Wrote this helper function to pipeline through a fastq file based on a position list and re-write it out to a compressed gz format
# This is a memory saving function BUT super slow (3-6X slower? than opening the file, doing vectorization operation and writing it)
# 161121: If we can do it in chunks, it's likely to work faster while still minimizing on memory draw. From testing with read.fastq.gz, 10^6 seems to work ok
write.fastq.gz = function(output.data, output.file, output.info) {
  
  lapply(c(paste0("\nWrite.fastq.gz started at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Input file: ", input.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Output file: ", output.file)), write, output.info, append=TRUE, ncolumns=1000)
  
  #input.con = file(input.file, "r")
  output.con = gzfile(output.file, "w")
  
  #line.list = unlist(scan(file=line.list.file, what=numeric()))
  
  curr.line = 1
  nlines = 1000000
  end.line = curr.line + nlines - 1
  
  # write data in chunks nlines 
  while(end.line <= length(output.data)) {
    
    cat(output.data[curr.line:end.line], file=output.con, sep="\n")
    
    curr.line = curr.line + nlines
    end.line = curr.line + nlines - 1
    
  }
  
  # There are two cases when we get here. Either there were exactly end.line lines OR there were not enough lines and we've skipped the last portion
  
  
  
  
  list.pos = 1
  end.list.pos = min((list.pos + nlines/4 - 1), length(line.list))
  lines.written = 0

  curr.line.list = unlist(lapply(line.list[list.pos:end.list.pos], function(x) seq(x, x+3, 1)))
  lapply(c(paste0("curr.line.list final number: ", as.character(tail(curr.line.list, 1)))), write, output.info, append=TRUE, ncolumns=1000)
  
  
  
  #while (length(curr.data <- scan(file=input.con, what=character(), nlines=4, sep="\n", quiet=TRUE)) > 0) {
  while (length(curr.data <- readLines(con=input.con, n=nlines)) > 3) {
    
    #lapply(c(paste0("curr.line: ", curr.line)), write, output.info, append=TRUE, ncolumns=1000)
    #lapply(c(paste0("lines read: ", length(curr.data))), write, output.info, append=TRUE, ncolumns=1000)
    
    # generate a line sequence for comparison against the list of read lines
    curr.data.lines = seq(curr.line, min((curr.line + nlines-1), (curr.line + length(curr.data) - 1)), 1)
    
    # grab any that match from the position set and write it to the output connection right away
    # Need to set up a conditional here in case there are no matches, otherwise it will output an extra \n
    if (length(curr.data[curr.data.lines %in% curr.line.list]) > 0) {
      cat(curr.data[curr.data.lines %in% curr.line.list], file=output.con, sep="\n")  
    }
    
    lines.written = lines.written + length(curr.data[curr.data.lines %in% curr.line.list])
    #lapply(c(paste0("Lines written: ", lines.written)), write, output.info, append=TRUE, ncolumns=1000)
    
    
    # Now check if we need to advance curr.line.list when our current position in the file is past the largest line number in our line list
    # if (tail(curr.data.lines, 1) >= tail(curr.line.list, 1)) {
    while (tail(curr.data.lines, 1) >= tail(curr.line.list, 1)) {
      
      # We need a special boundary case for when we are going past the line.list 
      if (end.list.pos >= length(line.list)) { break }
      
      #lapply(c(paste0("curr.data.lines final number: ", as.character(tail(curr.data.lines, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      #lapply(c(paste0("curr.line.list final number: ", as.character(tail(curr.line.list, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      
      # Boundary statement above should be okay. We should always be advancing this far otherwise
      list.pos = list.pos + nlines/4
      end.list.pos = min((list.pos + nlines/4 - 1), length(line.list))
      #lapply(c(paste0("new list.pos: ", as.character(list.pos))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      #lapply(c(paste0("end list.pos: ", as.character(end.list.pos))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      
      curr.line.list = unlist(lapply(line.list[list.pos:end.list.pos], function(x) seq(x, x+3, 1)))
      #lapply(c(paste0("new curr.line.list start number: ", as.character(head(curr.line.list, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      #lapply(c(paste0("new curr.line.list final number: ", as.character(tail(curr.line.list, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      
      # Now we need to recheck that segment for overlap that's outside the bounds and write any to the disk
      if (length(curr.data[curr.data.lines %in% curr.line.list]) > 0) {
        cat(curr.data[curr.data.lines %in% curr.line.list], file=output.con, sep="\n")  
      }
      lines.written = lines.written + length(curr.data[curr.data.lines %in% curr.line.list])
      #lapply(c(paste0("Lines written in loop: ", lines.written)), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
    }
    
#     if (curr.line == line.list[list.pos]) {
#       
#       cat(curr.data, file=output.con, sep="\n")
#       
#       
#       list.pos = list.pos + 1
#     }
    # advance to the next line set
    if (tail(curr.line.list, 1) < tail(curr.data.lines, 1)) { break }
    curr.line = curr.line + nlines
    
  }
  
  close(input.con)
  close(output.con)
  lapply(c(paste0("Lines written: ", lines.written)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Completed at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
}

# Wrote this helper function to pipeline through a fastq file based on a position list and re-write it out to a compressed gz format
# This is a memory saving function BUT super slow (3-6X slower? than opening the file, doing vectorization operation and writing it)
# 161121: If we can do it in chunks, it's likely to work faster while still minimizing on memory draw. From testing with read.fastq.gz, 10^6 seems to work ok
write.fastq.gz.old = function(input.file, output.file, line.list, output.info) {
  
  lapply(c(paste0("\nWrite.fastq.gz started at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Input file: ", input.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Output file: ", output.file)), write, output.info, append=TRUE, ncolumns=1000)
  
  input.con = file(input.file, "r")
  output.con = gzfile(output.file, "w")
  
  #line.list = unlist(scan(file=line.list.file, what=numeric()))
  
  curr.line = 1
  nlines = 1000000
  list.pos = 1
  end.list.pos = min((list.pos + nlines/4 - 1), length(line.list))
  lines.written = 0
  
  curr.line.list = unlist(lapply(line.list[list.pos:end.list.pos], function(x) seq(x, x+3, 1)))
  lapply(c(paste0("curr.line.list final number: ", as.character(tail(curr.line.list, 1)))), write, output.info, append=TRUE, ncolumns=1000)
  
  #curr.data.lines = seq(curr.line, nlines, 1)
  
  
  #while (length(curr.data <- scan(file=input.con, what=character(), nlines=4, sep="\n", quiet=TRUE)) > 0) {
  while (length(curr.data <- readLines(con=input.con, n=nlines)) > 3) {
    
    #lapply(c(paste0("curr.line: ", curr.line)), write, output.info, append=TRUE, ncolumns=1000)
    #lapply(c(paste0("lines read: ", length(curr.data))), write, output.info, append=TRUE, ncolumns=1000)
    
    # generate a line sequence for comparison against the list of read lines
    curr.data.lines = seq(curr.line, min((curr.line + nlines-1), (curr.line + length(curr.data) - 1)), 1)
    
    # grab any that match from the position set and write it to the output connection right away
    # Need to set up a conditional here in case there are no matches, otherwise it will output an extra \n
    if (length(curr.data[curr.data.lines %in% curr.line.list]) > 0) {
      cat(curr.data[curr.data.lines %in% curr.line.list], file=output.con, sep="\n")  
    }
    
    lines.written = lines.written + length(curr.data[curr.data.lines %in% curr.line.list])
    #lapply(c(paste0("Lines written: ", lines.written)), write, output.info, append=TRUE, ncolumns=1000)
    
    
    # Now check if we need to advance curr.line.list when our current position in the file is past the largest line number in our line list
    # if (tail(curr.data.lines, 1) >= tail(curr.line.list, 1)) {
    while (tail(curr.data.lines, 1) >= tail(curr.line.list, 1)) {
      
      # We need a special boundary case for when we are going past the line.list 
      if (end.list.pos >= length(line.list)) { break }
      
      #lapply(c(paste0("curr.data.lines final number: ", as.character(tail(curr.data.lines, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      #lapply(c(paste0("curr.line.list final number: ", as.character(tail(curr.line.list, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      
      # Boundary statement above should be okay. We should always be advancing this far otherwise
      list.pos = list.pos + nlines/4
      end.list.pos = min((list.pos + nlines/4 - 1), length(line.list))
      #lapply(c(paste0("new list.pos: ", as.character(list.pos))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      #lapply(c(paste0("end list.pos: ", as.character(end.list.pos))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      
      curr.line.list = unlist(lapply(line.list[list.pos:end.list.pos], function(x) seq(x, x+3, 1)))
      #lapply(c(paste0("new curr.line.list start number: ", as.character(head(curr.line.list, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      #lapply(c(paste0("new curr.line.list final number: ", as.character(tail(curr.line.list, 1)))), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
      
      # Now we need to recheck that segment for overlap that's outside the bounds and write any to the disk
      if (length(curr.data[curr.data.lines %in% curr.line.list]) > 0) {
        cat(curr.data[curr.data.lines %in% curr.line.list], file=output.con, sep="\n")  
      }
      lines.written = lines.written + length(curr.data[curr.data.lines %in% curr.line.list])
      #lapply(c(paste0("Lines written in loop: ", lines.written)), write, "161122.write.test.info", append=TRUE, ncolumns=1000)
    }
    
    #     if (curr.line == line.list[list.pos]) {
    #       
    #       cat(curr.data, file=output.con, sep="\n")
    #       
    #       
    #       list.pos = list.pos + 1
    #     }
    # advance to the next line set
    if (tail(curr.line.list, 1) < tail(curr.data.lines, 1)) { break }
    curr.line = curr.line + nlines
    
  }
  
  close(input.con)
  close(output.con)
  lapply(c(paste0("Lines written: ", lines.written)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Completed at: ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
}

# 150406: Similar to make.seq.files except takes into account that there are 4 files per samples as there are 4 lanes on the NextSeq
# 150403: Found out that the scan function can read in .gz files! No more pre-unzipping them. Just use as is.

make.NextSeq.files = function(MS.fastq.file, seq.directory, value.directory, start.row, end.row) {
  
  file.end = "_R1_001.fastq.gz"
  #file.end = "_R1_001.fastq"
  
  seq.directory = "seqFiles"
  value.directory = "seqValues"
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location
  #output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(directory, x), collapse="")))
  
  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], seq.directory), collapse="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], value.directory), collapse="/"))
  
  for (i in start.row:end.row) {
    
    for (j in 1:4) {
      
      # So you'll need to iterate through the series of 4 versions and save all the data into a single file at the end
      curr.fastq.file = paste0(c(MS.fastq.loc[i], j, file.end), collapse="")

      curr.fastq = scan(file=curr.fastq.file, what=(seq=""), sep="\n")
      curr.fastq.seq = curr.fastq[seq(2L, length(curr.fastq), 4L)] # grab lines starting at 2 and every 4th after that
      curr.fastq.values = curr.fastq[seq(4L, length(curr.fastq), 4L)]
      
      # Now be sure to tally the values instead of just replacing them into the table
      MS.fastq.data$seq.file.size[i] = as.numeric(MS.fastq.data$seq.file.size[i]) + length(curr.fastq.seq)
      MS.fastq.data$val.file.size[i] = as.numeric(MS.fastq.data$val.file.size[i]) + length(curr.fastq.values)

      # Write the sequence and value files based on the location listed in the original MS.fastq.file
      write.table(curr.fastq.seq, file=MS.fastq.data$seqLocation[i], append=TRUE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)      
      write.table(curr.fastq.values, file=MS.fastq.data$valLocation[i], append=TRUE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    }
  }
 
  # Write the values (reads per sample) out to the MS.fastq.file for comparison and checking later.
  write.table(MS.fastq.data, file=MS.fastq.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

}

### ----------------------- Demultiplex files. Could potentially update or change -----------------------###

# Demultiplex files using barcode sequences into files based on template names
# Techincally when Owen initially works with the fasta files, any failing clusters are not included BUT
# sequences with sporadic low values could still make it in through the process
# At this point we aren't investigating the quality value but it may be useful to look up later
# 141218: Need to fix this to be more efficient and line 2 = sequence, line 4 = values
demux.seq.files = function(MS.fastq.file, seq.fq.file, barcode.fq.file) {
  
  # Read both seq and barcode fastq.files
  # You can break the barcode file into just looking at the sequence data but you need to take all 4 lines of the seq.fq.data and put it into the new files
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location
  MS.fastq.bc = MS.fastq.data$barcode
   
  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1]), collapse="/"))
  
  seq.fq.data = scan(file=seq.fq.file, what=(seq=""), sep="\n")
  barcode.fq.data = scan(file=barcode.fq.file, what=(seq=""), sep="\n")

  # Reformat and pull barcodes off but do we do a single basepair mismatch algorithm? 
  # With a single mismatch algorithm, it generates 32 possible acceptable choices or 40 if we count "N" as acceptable
  
  for (i in 1:nrow(MS.fastq.data)) {
    
    curr.fq.file = MS.fastq.loc[i]
    curr.fq.bc = barcode.fq.data[i]
    curr.fq.bc = getComp(curr.fq.bc)
    
    # Enumerate all the possible barcodes
    all.curr.bc = enum.barcodes(curr.fq.bc)
    
    # Check the barcode file for matches - this will give all matching positions in the other file BUT it's the second line and every 4th thereafter
    # Use "which" to get the index value
    match.bc.lines = which(barcode.fq.data %in% all.curr.bc)
    
    # Now we need to extract the lines from X-1, X, X+1, X+2
    match.bc.line = match.bc.lines - 1
    
    # How do we extract all 4 lines? make a sequence vector with lapply
    match.lines = unlist(lapply(match.bc.line, function(x) seq(x, x+3, 1)))
    
    # Now just grab all the appropriate lines
    curr.fastq.data = seq.fq.data[match.lines[1]]
    write.table(curr.fastq.data, file=curr.fq.file, append=FALSE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    for (j in 2:length(match.lines)) {
      
      write.table(seq.fq.data[match.lines[j]], file=curr.fq.file, append=TRUE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
      
    }
    
  }

}


demux.seq.files.v2 = function(MS.fastq.file, seq.fq.file, barcode.fq.file) {
  
  # Read both seq and barcode fastq.files
  # You can break the barcode file into just looking at the sequence data but you need to take all 4 lines of the seq.fq.data and put it into the new files
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location
  MS.fastq.bc = MS.fastq.data$barcode
  
  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1]), collapse="/"))
  
  seq.fq.data = scan(file=seq.fq.file, what=(seq=""), sep="\n")
  barcode.fq.data = scan(file=barcode.fq.file, what=(seq=""), sep="\n")
  # grab every second line - barcode sequences only
  barcode.fq.data = barcode.fq.data[seq(2L, length(barcode.fq.data), 4L)]
  
  # Reformat and pull barcodes off but do we do a single basepair mismatch algorithm? 
  # With a single mismatch algorithm, it generates 32 possible acceptable choices or 40 if we count "N" as acceptable
  # first version was very slow and might be running out of memory. Need to be slower code but more memory efficient
  
  # This barcode.list will store the enumerated possibilities of each barcode
  barcode.list = list()
  
  # enumerate all versions of the barcodes
  for (i in 1:length(MS.fastq.bc)) {
    
    barcode.list[i] = enum.barcodes(MS.fastq.bc[i])    
  }
  
  
  for (i in 1:length(barcode.fq.data)) {
    
    start.line = (i-1)*4+1
    
    # grab the current 4 lines
    curr.entry = seq.fq.data[seq(start.line, start.line+4, 1)]
    
    # grab the corresponding barcode
    curr.barcode = barcode.fq.data[i]

    for (j in 1:length(barcode.list)) {
      
      if (curr.barcode %in% barcode.list[[j]]) {
        
        write.table(curr.entry, file=MS.fastq.loc[j], append=TRUE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
        break
      }
      
    }
      
  }
  
}

# Read the data in line by line
demux.seq.files.v3 = function(MS.fastq.file, seq.fq.file, barcode.fq.file) {
  
  # Read both seq and barcode fastq.files
  # You can break the barcode file into just looking at the sequence data but you need to take all 4 lines of the seq.fq.data and put it into the new files
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location
  MS.fastq.bc = MS.fastq.data$barcode
  
  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1]), collapse="/"))
  
  #seq.fq.data = scan(file=seq.fq.file, what=(seq=""), sep="\n")
  #barcode.fq.data = scan(file=barcode.fq.file, what=(seq=""), sep="\n")
  # grab every second line - barcode sequences only
  #barcode.fq.data = barcode.fq.data[seq(2L, length(barcode.fq.data), 4L)]
  
  # Reformat and pull barcodes off but do we do a single basepair mismatch algorithm? 
  # With a single mismatch algorithm, it generates 32 possible acceptable choices or 40 if we count "N" as acceptable
  # first version was very slow and might be running out of memory. Need to be slower code but more memory efficient
  
  # This barcode.list will store the enumerated possibilities of each barcode
  barcode.list = list()
  
  # enumerate all versions of the barcodes
  for (i in 1:length(MS.fastq.bc)) {
    
    barcode.list[i] = list(enum.barcodes(MS.fastq.bc[i]))    
  }
  
  seq.conn = file(seq.fq.file, open="r")
  bar.conn = file(barcode.fq.file, open="r")
  
  while(length(curr.entry <- readLines(seq.conn, n=4, ok=TRUE)) == 4) {
    
    # curr.entry = readLines(seq.conn, n=4, ok=TRUE)
    curr.barcode = readLines(bar.conn, n=4, ok=TRUE)
    if (length(curr.barcode) < 4) { break }
    
    # grab the corresponding barcode
    curr.barcode = curr.barcode[2]
    
    for (j in 1:length(barcode.list)) {
      
      if (curr.barcode %in% barcode.list[[j]]) {
        
        write.table(curr.entry, file=MS.fastq.loc[j], append=TRUE, sep="\n", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
        break
      }
      
    }
    
  }
  
  close(seq.conn)
  close(bar.conn)
  
}

### ----------------------- Helper Functions not directly part of the initial analysis -----------------------###

# This is a quick analysis tool to help in identifying the spread/population of gap-reads that exist after the ligation arm read. 
# The final number does not likely take into account repeats of the same barcode although this is currently a minimal amount (5% or so)
MIP.gap.summary = function(MS.fastq.file, MIP.info.file, output.directory, MIP.pos.list, start.row, end.row) {

  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  
  # Create the output directory
  dir.create(output.directory)
  
  # create a table of all the fastq file names
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.data = MS.fastq.data[start.row:end.row,]
  
  # generate a list of all the fastq file locations
  MS.fastq.loc = MS.fastq.data$seqLocation

  # experiment names for file names later
  MS.exp.name = MS.fastq.data$exp.name

  # MIP information table - the same used for other analyses
  MIP.info = read.table(MIP.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  for (i in 1:length(MS.fastq.loc)) {
    
    # Read in the given fastq sequence data file 
    curr.fastq.seq = as.list(scan(file=MS.fastq.loc[i], what=(seq="")))
    
    for (j in 1:length(MIP.pos.list)) {
      
      curr.lig.arm = toupper(MIP.info$lig.seq.read[MIP.pos.list[j]])
      curr.MIP = MIP.info$MIP.name[MIP.pos.list[j]]
      #curr.lig.arm = substr(curr.lig.arm, 1, (nchar(curr.lig.arm)-5))
            
      # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
      #matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, curr.fastq.seq)]
      
      # 151209: Modified to only count and show data from exact ligation arm matches in the correct position
      matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, lapply(curr.fastq.seq, function(x) substr(x, MIP.barcode.length+1, MIP.barcode.length+nchar(curr.lig.arm))), perl=TRUE)]
      
      # Convert it to a table with frequencies but make it a data frame first
      matching.lig.seq = as.data.frame(unlist(matching.lig.seq))
      
      # Making the table in this fashion will only leave frequencies but not by unique bacodes
      #matching.gap = as.data.frame(table(lapply(matching.lig.seq, function(x) substr(x, 12+nchar(curr.lig.arm)+1, 50))))
      matching.gap = as.data.frame(table(lapply(matching.lig.seq, function(x) substr(x, 12+nchar(curr.lig.arm)+1, 75))))
      
      matching.gap = matching.gap[order(-matching.gap[2], matching.gap[1]),]
      
      # Save the file
      outputFile = paste(c(output.directory, paste(c(MS.exp.name[i], curr.MIP, curr.lig.arm, "gap.fill.seq.freq.txt"), collapse = ".")), collapse="")
      outputFile = paste(c(output.directory, paste(c(MS.exp.name[i], curr.MIP, "gap.fill.seq.freq.txt"), collapse = ".")), collapse="")
      
      write.table(matching.gap, file=outputFile, append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
      
    }
  }
}

# This is a quick analysis tool to help in identifying the spread/population of gap-reads that exist after the ligation arm read. 
# The final number does not likely take into account repeats of the same barcode although this is currently a minimal amount (5% or so)
# 160223: Generated a new version because it looks like the rDNA locus may house a high % (~10)
# Generate an additional column on if the 2 changes found place the gap sequence as wt vs mut (plasmid)
# Need to only grab samples with a
rDNA.MIP.gap.summary = function(MS.fastq.file, MIP.info.file, output.directory, MIP.pos.list, start.row, end.row) {
  
  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  
  # Create the output directory
  dir.create(output.directory)
  
  # create a table of all the fastq file names
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.data = MS.fastq.data[start.row:end.row,]
  
  # generate a list of all the fastq file locations
  MS.fastq.loc = MS.fastq.data$seqLocation
  MS.val.loc = MS.fastq.data$valLocation
  
  # experiment names for file names later
  MS.exp.name = MS.fastq.data$exp.name
  
  # MIP information table - the same used for other analyses
  MIP.info = read.table(MIP.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  for (i in 1:length(MS.fastq.loc)) {
    
    # Read in the given fastq sequence data file 
    curr.fastq.seq = as.list(scan(file=MS.fastq.loc[i], what=(seq="")))
    curr.fastq.val = as.list(scan(file=MS.val.loc[i], what=(seq="")))
    for (j in 1:length(MIP.pos.list)) {
      
      curr.lig.arm = toupper(MIP.info$lig.seq.read[MIP.pos.list[j]])
      curr.MIP = MIP.info$MIP.name[MIP.pos.list[j]]
      curr.MIP.lig.gap = as.numeric(unlist(strsplit(MIP.info$lig.gap[MIP.pos.list[j]], split=";")))
      
      # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
      #matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, curr.fastq.seq)]
      
      # 151209: Modified to only count and show data from exact ligation arm matches in the correct position
      matching.lig.entries = grep(curr.lig.arm, lapply(curr.fastq.seq, function(x) substr(x, MIP.barcode.length+1, MIP.barcode.length+nchar(curr.lig.arm))))
      
      #matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, lapply(curr.fastq.seq, function(x) substr(x, MIP.barcode.length+1, MIP.barcode.length+nchar(curr.lig.arm))))]
      #matching.lig.seq = as.data.frame(unlist(curr.fastq.seq[matching.lig.entries]))
      #matching.lig.val = as.data.frame(unlist(curr.fastq.val[matching.lig.entries]))
      
      matching.lig.seq = (unlist(curr.fastq.seq[matching.lig.entries]))
      matching.lig.val = (unlist(curr.fastq.val[matching.lig.entries]))
      
      # Grab just the gap-fill sequences
      matching.gap.seq = unlist(lapply(matching.lig.seq, function(x) substr(x, 12+nchar(curr.lig.arm)+1, 50)))
      matching.gap.val = unlist(lapply(matching.lig.val, function(x) substr(x, 12+nchar(curr.lig.arm)+1, 50)))
      
      # Now remove any entries that are less than min.quality along the gap.fill sequence
      # matching.gap.seq = matching.gap.seq[!unlist(lapply(matching.gap.val, function(x) any(unlist(get.quality.score(x)) < 20)))]
      
      # Convert it to a table with frequencies but make it a data frame first
      matching.gap = as.data.frame(table(matching.gap.seq))
      for (k in 1:length(curr.MIP.lig.gap)) {
        
        matching.gap[,(ncol(matching.gap)+1)] = unlist(lapply(matching.gap$matching.gap.seq, function(x) substr(x, (curr.MIP.lig.gap[k]+1), (curr.MIP.lig.gap[k]+1))))
      }
      
      
      # Making the table in this fashion will only leave frequencies but not by unique bacodes
      #matching.gap = as.data.frame(table(lapply(matching.lig.seq, function(x) substr(x, 12+nchar(curr.lig.arm)+1, 50))))
      
      matching.gap = matching.gap[order(-matching.gap[2], matching.gap[1]),]
      
      # Save the file
      #outputFile = paste(c(output.directory, paste(c(MS.exp.name[i], curr.MIP, curr.lig.arm, "gap.fill.seq.freq.txt"), collapse = ".")), collapse="")
      outputFile = paste(c(output.directory, paste(c(MS.exp.name[i], curr.MIP, "gap.fill.seq.freq.txt"), collapse = ".")), collapse="")
      
      write.table(matching.gap, file=outputFile, append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
      
    }
  }
}


# An original version of MIP.gap.summary
make.MIP.table = function(MIP.info.file, MS.fastq.file, MIP.nums, outputFile){
  
  # Open up the MIP.info table
  MIP.info = read.table(MIP.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  # scan in the specific file
  curr.fastq.seq = as.list(scan(file=MS.fastq.file, what=(seq="")))
  
  for (i in 1:length(MIP.nums)){
    
    curr.lig.arm = toupper(MIP.info$lig.seq.read[MIP.nums[i]])
    curr.read.wt = toupper(MIP.info$wt.gap.fill.read[MIP.nums[i]])
    curr.read.snv = toupper(MIP.info$snv.gap.fill.read[MIP.nums[i]])
    
    # WT lig + gap.fill read
    curr.MIP.wt = paste(c(curr.lig.arm, curr.read.wt), collapse="")
    # SNV lig + gap.fill read
    curr.MIP.snv = paste(c(curr.lig.arm, curr.read.snv), collapse="")
    
    # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
    matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, curr.fastq.seq, perl=TRUE)]
    
    # Subset into barcodes and reads
    #matching.lig.barcodes = unlist(lapply(matching.lig.seq, function (x) substr(x, 1, 12)))
    #matching.lig.reads = unlist(lapply(matching.lig.seq, function(x) substr(x, 13,50)))
    
    gap.fill.start = 12+nchar(curr.lig.arm)+1
    matching.lig.gaps = unlist(lapply(matching.lig.seq, function(x) substr(x, gap.fill.start, MIP.read.length)))
    #   matching.lig.gaps = unlist(lapply(matching.lig.seq, function(x) substr(x, 12, MIP.read.length)))
    
    #matching.lig.gaps.table = table(sort(unique(matching.lig.gaps)))
    matching.lig.gaps.table = as.data.frame(table(matching.lig.gaps))
    matching.lig.gaps.table = matching.lig.gaps.table[order(-matching.lig.gaps.table$Freq),]
    write.table(matching.lig.gaps.table, file=paste(c(MIP.nums[i],outputFile), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)  
    
    return (matching.lig.gaps.table)
    
  }
}

# A function to help generate MIP.info files for different pools
# Takes in a list of all "working" MIP data and strain data. If the working.status = yes then it will keep those in the final info file

make.MIP.info.table = function(all.MIP.info.file, strain.file, output.file) {
  
  all.MIP.info = read.table(all.MIP.info.file, header=TRUE, sep="\t", colClasses="character")
  strain.list = read.table(strain.file, header=TRUE, sep="\t", colClasses="character")
  
  new.MIP.info = all.MIP.info[all.MIP.info$snv.strain %in% strain.list$strain,]
  new.MIP.info = new.MIP.info[new.MIP.info$working.status == "yes",]
  
  write.table(new.MIP.info, file=output.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  return (new.MIP.info)
  
}

### ----------------------- Derelict functions that appear to no longer have function -----------------------###

# MS.run is the current location of the file we are interested in analysing
getMatchingGaps = function(MS.run.file, MIP.info) {

  curr.MIP.list = list()
  
  
  # Generate the list of ligation arms that you want analysed
  MIP.lig.arms = unlist(lapply(MIP.info$ligation.arm.read, toupper))  
  
  for (i in 1:nrow(MIP.info)) {
    
#     curr.MIP.file = paste(c(MS.run, MIP.lig.arms[i], "txt"), collapse=".")
#     
#     curr.MIP = read.table(curr.MIP.file, header=FALSE, sep="", colClasses = "character")
#     
#     MIP.list$MiSeq.130329[i] = sum(as.numeric(curr.MIP[[1]]))

    
    
    if (sum(curr.MIP[,2] %in% toupper(MIP.list$wt.gap.fill.read[i])) == 1) {  
      MIP.list$WT.total[i] = curr.MIP[,1][curr.MIP[,2] %in% toupper(MIP.list$wt.gap.fill.read[i])]
    }
    
    if (sum(curr.MIP[,2] %in% toupper(MIP.list$SNV.gap.fill.read[i])) == 1) {  
      MIP.list$SNV.total[i] = curr.MIP[,1][curr.MIP[,2] %in% toupper(MIP.list$SNV.gap.fill.read[i])]
    }
    
    curr.MIP.list = c(curr.MIP.list, curr.MIP)        
  }
  
  return (MIP.list)
  

}


sumMatchingReads = function(MS.run, MIP.lig.arms, MIP.list) {
  
  total.matching.reads = 0
  
  for (i in 1:nrow(MIP.list)) {
    
    curr.MIP.file = paste(c(MS.run, MIP.lig.arms[i], "txt"), collapse=".")
    
    curr.MIP = read.table(curr.MIP.file, header=FALSE, sep="", colClasses = "character")
    
    total.matching.reads = total.matching.reads + sum(as.numeric(curr.MIP[,1]))
    
  }
  
  return (total.matching.reads)  
  
}



analyse.MIPS = function(MS.fastq.file, MIP.info.file, directory) {

  #results = analyse.MIPS(MS.fastq.file, MIP.info.file, "./130411_MiSeq/analysis2/")
#   MS.fastq.file = "130329.MS.fastq.files.txt"
#   MIP.info.file = "MS130329.all.counts.txt"
#   directory = "./test1/"
  dir.create(directory)
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")  
  MS.fastq.loc = MS.fastq.data$seqLocation
  output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(directory, x), collapse="")))
  
  MIP.info = read.table(MIP.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  MIP.data = list()
  
  for (i in 1:length(MS.fastq.loc)) {
    
    MIP.data[[i]] = get.MIP.Barcodes(MIP.info, MS.fastq.loc[i])    
    # Now save this table to the running directory      
    write.table(MIP.data[[i]], file=paste(c(output.prefix.list[i], "results", "txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    
  }
  
  
  # Now collate all the data into a larger file of information
  all.MIP.info = summarize.MIPS(MIP.info, MIP.data)
  write.table(all.MIP.info, file=paste(c(directory, "MIP.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  
  return (list(MIP.data, all.MIP.info))
}

# -------------------- Re-analysis of MIPS scripts -------------------- #

# 150326: The purpose of this series of scripts is to re-analyse data from expName.results.txt in multiple files.
# Given a secondary file of strains to remove, it would then proceed to find the strain.MIPs, and calculate an overall abundance for average.
# From there, it would have to summ that abundance across all the strains to remove.
# This % would then be converted to an expected reads number for each remaing MIP based on WT+SNV.
# This expected number would be removed from the WT counts of that MIP and a SNV/Total abundance would be recalculated.
# From there we can recalculate the average abundance of the remaining strains

recalculate.abundance = function(MIP.avg.file, output.file, rm.strains.file, output.info) {
  
  lapply(c(paste0("Recalculate abundance called at ", Sys.time())), write, output.info, append=FALSE, ncolumns=1000)
  lapply(c(paste0("MIP.avg.file: ", MIP.avg.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("output.file: ", output.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("rm.strains.file: ", rm.strains.file)), write, output.info, append=TRUE, ncolumns=1000)
  
  rm.strains = read.table(rm.strains.file, header=TRUE, sep="\t", colClasses = "character")$strain  
  
  lapply(c(paste0("strains to remove: ", rm.strains)), write, output.info, append=TRUE, ncolumns=1000)
  
  MIP.avg.df = read.table(MIP.avg.file, row.names=1, header=TRUE, sep="\t")
 
  # Grab the relevant strain data into a remove and keep set
  rm.avg.df = MIP.avg.df[rownames(MIP.avg.df) %in% rm.strains,]
  keep.avg.df = MIP.avg.df[!rownames(MIP.avg.df) %in% rm.strains,]
  
  # Sum the average abundance of rm.strains
  rm.avg.df.sum = colSums(rm.avg.df)
  # Generate a new denominator to divide into the values
  new.avg.denom = 1-rm.avg.df.sum
  
  # Now divide it into each column value with apply
  new.avg.df = t(apply(keep.avg.df, 1, function(x) x/new.avg.denom))
  
  write.table(new.avg.df, file=output.file, append=FALSE, sep="\t", eol = "\n", row.names=TRUE, col.names=TRUE, quote=FALSE)  
  
}



# recalculate.abundance = function(input.dir, output.dir, MS.fastq.file, remove.strains.file){
#   
#   input.dir = paste0(c("./", input.dir, "/"), collapse="")
#   output.dir = paste(c("./", output.directory, "/"), collapse = "")
#   output.info = paste(c(output.dir, "run.info.txt"), collapse = "")
#   
#   lapply(c(paste0("Recalculate abundance called at", Sys.time())), write, output.info, append=FALSE, ncolumns=1000)
#   lapply(c(paste0("input.dir: ", input.dir)), write, output.info, append=TRUE, ncolumns=1000)
#   lapply(c(paste0("output.dir: ", output.dir)), write, output.info, append=TRUE, ncolumns=1000)
#   lapply(c(paste0("MS.fastq.file: ", MS.fastq.file)), write, output.info, append=TRUE, ncolumns=1000)
#   lapply(c(paste0("remove.strains.file: ", remove.strains.file)), write, output.info, append=TRUE, ncolumns=1000)
#   
#   lapply(c(paste0("Total experiments: ", nrow(MS.fastq.data))), write, output.info, append=TRUE, ncolumns=1000)
#   lapply(c(paste0("Total MIPS to evaluate: ", nrow(MIP.info))), write, output.info, append=TRUE, ncolumns=1000)
#   
#   
#   
#   
#   file.names = read.table(exp.name.file, header=TRUE, sep="\t", colClasses = "character")$exp.name
#   remove.strains = read.table(remove.strains.file, header=TRUE, sep="\t", colClasses = "character")$strain  
#   
#   for (i in 1:length(file.names)) {
#     
#     # Open the file
#     curr.input.file = paste0(c(input.dir, file.names$exp.name, ".results.txt"), collapse="")
#     curr.input = read.table(curr.input.file, header=TRUE, sep="\t", colClasses = c(rep("character", 39), rep("numeric", 2), rep("character", 14)))
#     
#     # Relevant columns: $snv.strain, $wt.bc.unique, $snv.bc.unique, $snv.wt.unique.ratio
#     # Find the strains to remove
#     
#     lapply(c(paste0("Current file: ", curr.input.file)), write, output.info, append=TRUE, ncolumns=1000)
#     
#     removal.percent = 0
#     
#     # Go through each strain and recalculate it's abundance in that file
#     for (j in 1:length(remove.strains)) {      
#     }
#   }  
# }