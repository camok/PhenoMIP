# This is a supplementary script used in performing correlations 
# and generating heatmaps of the multiplex pooling data
# 150511: Version 6 updated to use "Generation" information instead of "Expansions"
# Basic assumption: E01 = Generation 1-3. Then G03 is moved to E02 which end at G04. 
# Nomenclature uses LAST generation on plate and assumes it will be more represented (L1s vs Adults)
# 151012: Version 7 updated to do an additional normalization by a strain's "average" curve. 
# This way the curve graphs will perhaps be tighter OR it may be easier to identify deviation?
# Also want to update with additional graphs to look at deviations at specific time points like a kind of box and whisker plot?
# Maybe re-instate the correlation plots to find replicates that behave similarly.
# 180426: Update analysis to focus on fold-change between each generation. This should generate an average fold change per strain per condition
# 180718: v10 update to change the minimum threshold for existence/abundance to 0.000541 (calculated mean of 0.000127 + 3*0.000138) will capture 99.7% of possible false positives


# Correlogram
library(lattice)
library(reshape2)
library(ggplot2)
library(scales)
library(directlabels)
library(RColorBrewer)
library(statmod)
library(EnvStats)
library(ggrepel)
library(gplots)

#options(java.home="C:\\Program Files\\Java\\jdk1.8.0_201/")
Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk-11.0.2\\")
library(xlsx)

## ---------- Global variables ---------- ##
font.size=24

# Note for "range" in a boxplot, default = 1.5IQR, strong outliers = 3IQR, 0 = whiskers at extreme datapoints (no outliers)
# The whiskers defined by the last data point within the range, all others points are outliers. This is why it can look uneven.
IQ.range = 1.5

output.info = "run.info.txt"
combined.data.file = "combined.data.xlsx"

# Assume E01 ends at G03. E02 = G04, etc...
# NEW: E01 will be G0 as it was not on RNAi, then the first expansion on RNAi (E02) = G01
gen.offset = -1

# Some assumptions about the data:
# What is the minimum value we can expect without being 0? 
# Given that 0 is a completely realistic value, if we are doing log transformations, can we convert it to something else for our purposes?
#min.abundance.value = 0.000542
#min.abundance.value = 0.000821
min.abundance.value = 0.001
# This number can depend a lot on whether or not we are log transforming data.
min.avg.rate = log2(min.abundance.value)

# What is the minimum start value we want to work with for population abundance in Generation 1? Is it worth analysing anything below 0.0025?
min.start.abundance.value = 0.0025

# Global variables to make passing data upwards easier
fc.rate.melt = list()
fc.total.melt = list()

## ---------- End Global variables ---------- ##

# 180505: Amalgamated function for completing multiple analyses
analyse.multiple.pools = function(dataset.file, output.dir, input.dir, row.start = 1, row.end=NULL) {
  

  dataset.df = read.table(dataset.file, header=TRUE, sep="\t", colClasses = c("character", "character", "character", "numeric", "numeric", "character", "character", "numeric"))
  run.command = paste0("analyse.multiple.pools(\"",dataset.file, "\", \"", output.dir, "\", \"", input.dir, "\", row.start=", row.start, ", row.end=", row.end, ")")
  
  #fold.change.zscore.list = list()
  #fold.change.zscore.pval.list = list()
  dir.create(output.dir, recursive=TRUE)
  run.info = paste0(c(output.dir,"/run.info.txt"), collapse="")
  
  lapply(c(paste0("MIP.graph.analysis.v9.R")), write, run.info, append=TRUE, ncolumns=1000)
  lapply(c(run.command), write, run.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Minimum starting abundance: ", min.start.abundance.value)), write, run.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Minimum abundance value: ", min.abundance.value)), write, run.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Program start at ", Sys.time())), write, run.info, append=TRUE, ncolumns=1000)
  
  if (is.null(row.end)) { row.end = nrow(dataset.df)}
  
  for (i in row.start:row.end) {
  
    input.info = paste0(c(input.dir, "/", dataset.df$data.file[i]), collapse="")
    output.data = paste0(c(output.dir, "/", dataset.df$output.dir[i], "/"), collapse="")
    group.info = paste0(c(input.dir, "/", dataset.df$group.file[i]), collapse="")
   
    lapply(c(paste0(c("Running analysis on: ", input.info), collapse="")), write, run.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0(c("Outputting to directory: ", output.data), collapse="")), write, run.info, append=TRUE, ncolumns=1000)

    generate.plots(input.info, output.data, dataset.df$num.rep[i], dataset.df$num.exp[i], group.info, dataset.df$norm.strain[i])

    
  }
  
}


# 141221: New function to make plots for our data. Instead of having to tailor the input each time, we should just 
# tell the program what the data initially looks like, and have it organize for the output for us.
# Information we need: data.file, num.replicates, num.expansions, group.file
# We'll have the data produce graphs based on grouping by expansion and replicates, but also colour the graphs by the group.file
# data.file: $strain, $XX.R0.E01, $XX.R1.E02, $XX.R2.E02, ..., $XX.Rn.Em ** 
# R00.E01 does not count as an expansion for num.exp and we only count sampled expansions, not total expansions
# group.file: $rep.num, $rep.name, $group 

generate.plots = function(data.file, output.directory, num.rep, num.exp, group.file=NULL, norm.strain="VC20019") {
  
  run.command = paste0(c("generate.plots(\"", data.file, "\", \"", output.directory, "\", ", num.rep, ", ", num.exp, ", \"", group.file, "\", \"", norm.strain, "\""), collapse="")
  
  output.directory = paste0(c("./", output.directory, "/"), collapse="")
  # Create the output directory
  dir.create(output.directory, recursive=TRUE)
  
  output.info <<- paste0(c(output.directory, "run.info.txt"), collapse = "")
  
  lapply(c(paste0("MIP.graph.analysis.v9.R")), write, output.info, append=FALSE, ncolumns=1000)
  lapply(c(run.command), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Program start at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  # When you read the table in, if you don't specify a colClasses type, it looks like it will guess at the format. 
  # Fortunately it looks like it takes our data and makes it numeric. The row.names call makes the first column into our row.names and removes it from the table
  data.df.exp = read.table(data.file, header=TRUE, sep="\t", row.names=1)
  #data.norm.df = as.data.frame(normalize.strains(data.df, norm.strain))
 
  # ---------- Data TRIM ---------- #
  # 180716: An issue that has been cropping up in analysis is that of the initial data in earlier poolings. 
  # What happens the first expansion has data for a sample that is below the false positive threshold
  # Data with sample of that format should be trimmed away and removed from further analysis
  # 1) Identify strains below threshold
  # 2) Verify that the strain has no real representation across the population
  # 3 Remove that entire row from the data frame
  
  trimmed.strains = which(data.df.exp[,1] < min.start.abundance.value)
  lapply(c(paste0("Total strains below start threshold: ", length(trimmed.strains))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0(rownames(data.df.exp)[trimmed.strains])), write, output.info, append=TRUE, ncolumns=1000)

  if (length(trimmed.strains) > 0) {
    data.df.exp <- data.df.exp[-trimmed.strains,]
    lapply(c(paste0("Total strains remaining for analysis: ", nrow(data.df.exp))), write, output.info, append=TRUE, ncolumns=1000)
    # AT this point we will not re-adjust the abundance since the presumption is that the strains were not there anyways
  }
  
  # 180525: If we know the false positive rate of the MIPs we can generate a cutoff value for the purposes of our downstream analysis
  # was 0.0122% +/- 0.0129 from Mok et al. 2017. 
  # Latest paper suggests we can reliably detect at 0.000541 and above (3 standard deviations above the mean FP rate for our combined data sets)
  data.df.exp[data.df.exp < min.abundance.value] <- min.abundance.value
  
  # What do we need at this point? A list of Replicate names, and a list of expansion names
  # We can glean it from the column names in data.df
  col.info = unlist(strsplit(colnames(data.df.exp), split="[.]"))
  
  exp.names = unique(as.numeric(gsub('[A-Z]', "",unique(col.info[grep("E", col.info)])))) + gen.offset
  # rep.names = unique(col.info[grep("R", col.info)])
  rep.names = unique(as.numeric(gsub('[A-Z]', "",unique(col.info[grep("R", col.info)]))))
  lapply(c(paste0("Expansion names ", exp.names)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Replicate names ", rep.names)), write, output.info, append=TRUE, ncolumns=1000)

  # Check the group information. Is it null or are we given a file?
  if (class(group.file) == "NULL") {
    # We should generate a default type of group.data where each replicate is numbered 1-n, the groups are all "replicate" or "avg"
    group.data = data.frame(rep.num =as.character(seq(1,length(set.col)-1)), rep.name = as.character(seq(1,length(set.col)-1)), group = rep("replicate", length(set.col)-1), stringsAsFactors=FALSE)
  } else {
    group.data = read.table(group.file, header=TRUE, sep="\t")   
  }
  

  # ---------- First pass analysis will be done on the expansions as groups ---------- #
  # The data will initially come grouped by expansion (given our definition) and each expansion should have the same number of replicates
  
  #exp.col = seq(2,ncol(data.df.exp), by=floor(ncol(data.df.exp)/num.rep))
  exp.col = seq(2,ncol(data.df.exp), by=num.rep)
  
  # Is this correct? Number of expansions is incorrect...
  # make.exp.plot.data(data.df.exp, paste0(c(output.directory, "by.gen"), collapse=""), norm.strain, seq(2,ncol(data.df.exp), by=floor(ncol(data.df.exp)/num.rep)), exp.names, group.data, num.rep)
  make.exp.plot.data(data.df.exp, paste0(c(output.directory, "by.gen"), collapse=""), norm.strain, exp.col, exp.names, group.data, num.rep)
  lapply(c(paste0("by.exp completed at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  # ---------- Second pass analysis will be done on the replicates as groups ---------- #
  
  # First grab the row.names
  data.df.rep = data.df.exp[,0]
  data.df.rep = cbind(data.df.rep, data.df.exp[,1])
  data.df.rep.names = colnames(data.df.exp)[1]
  
  for (i in 2:(num.rep+1)) {
    
    # Based on num.rep, grab every num.rep column at 2,3,4, etc.
    curr.col = seq(i, ncol(data.df.exp), by = num.rep)
    data.df.rep = cbind(data.df.rep, data.df.exp[,curr.col])
    data.df.rep.names = c(data.df.rep.names, colnames(data.df.exp)[curr.col])
    
  }
  
  colnames(data.df.rep) = data.df.rep.names
  
  # First attempt at generating replicate data
  all.rep.results = make.rep.plot.data(data.df.rep, paste0(c(output.directory, "by.rep"), collapse=""), norm.strain, num.rep, num.exp, exp.names, group.data)

  lapply(c(paste0("by.rep completed at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  # it would be nice to break things into groups to generate the data if needed.
  # Based on initial data, it looks like groups are ordered ALPHABETICALLY when plotted, but the decision is not transparent
  # Take the extra step here to sort the group names
  if (class(group.file) != "NULL") {
    
    group.names = as.character(sort(unlist(unique(group.data$group))))
        
    #group.rep.results = list()
    
    #group.avg.df = subset(all.rep.results[[1]][[1]], select=c(1,2,(ncol(all.rep.results[[1]][[1]])-num.exp):ncol(all.rep.results[[1]][[1]])))
    #group.norm.avg.df = subset(all.rep.results[[2]][[1]], select=c(1,2,(ncol(all.rep.results[[2]][[1]])-num.exp):ncol(all.rep.results[[2]][[1]])))
    
    group.avg.df = subset(all.rep.results[[1]][[1]], select=c(1,2))
    group.norm.avg.df = subset(all.rep.results[[2]][[1]], select=c(1,2))
    
    avg.group.names = list()
    # We'll dissect data.df.rep into the subgroups and run the same analyses on them
    for (i in 1:length(group.names)) {
            
      curr.group = group.names[i]
      lapply(c(paste0("by.rep.", as.character(curr.group), " started at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
      
      # This is the "rep" number we should be using when pulling data from the data.df.rep
      curr.group.rows = which(group.data$group == curr.group)
      
      if (length(curr.group.rows) > 1){
      
        avg.group.names = c(avg.group.names, curr.group)
        curr.group.df = group.data[group.data$group == curr.group,]
        
        
        group.data.df.rep = subset(data.df.rep, select=c(1))
        
        # Now we need to re-assemble the data.df.rep by the replicates we are truly interested in. 
        # 150327: If gaps exist in our data set (like a rep has been removed) then the rep number won't necessarily match it's position
        # Instead we need to base it on the row number of the rep!
        for (j in 1:nrow(curr.group.df)) {
          
          # take the current replicate
          curr.rep = curr.group.rows[j]
          curr.rep.start = 2 + (curr.rep-1)*num.exp
          curr.rep.end = curr.rep.start + num.exp - 1
          
          lapply(c(paste0("group col selection is ", curr.rep.start, " to ", curr.rep.end)), write, output.info, append=TRUE, ncolumns=1000)
          
          # Subset allows us to keep the columns named properly
          group.data.df.rep = cbind(group.data.df.rep, subset(data.df.rep, select=c(curr.rep.start:curr.rep.end)))
          
        }
        
        # Each returned result will be made of the same strains. The removal of "0" presence strains at the beginning will be uniform across all groups
        # since this information is gleaned from the start population
        #group.rep.results[i] = make.rep.plot.data(data.df.rep, paste0(c(output.directory, "by.rep.", as.character(curr.group)), collapse=""), norm.strain, nrow(curr.group.df), num.exp, exp.names, curr.group.df)
        
        group.rep.results = make.rep.plot.data(group.data.df.rep, paste0(c(output.directory, "by.rep.", curr.group), collapse=""), norm.strain, nrow(curr.group.df), num.exp, exp.names, curr.group.df)
        
        group.avg.df = cbind(group.avg.df, subset(group.rep.results[[1]][[1]], select=c((ncol(group.rep.results[[1]][[1]])-num.exp):ncol(group.rep.results[[1]][[1]]))))
        group.norm.avg.df = cbind(group.norm.avg.df, subset(group.rep.results[[2]][[1]], select=c((ncol(group.rep.results[[2]][[1]])-num.exp):ncol(group.rep.results[[2]][[1]]))))
        
        lapply(c(paste0("by.rep.", as.character(curr.group), " completed at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
        
      }
    }
    
    group.avg.df = cbind(group.avg.df, subset(all.rep.results[[1]][[1]], select=c((ncol(all.rep.results[[1]][[1]])-num.exp):ncol(all.rep.results[[1]][[1]]))))
    group.norm.avg.df = cbind(group.norm.avg.df, subset(all.rep.results[[2]][[1]], select=c((ncol(all.rep.results[[2]][[1]])-num.exp):ncol(all.rep.results[[2]][[1]]))))


    group.names = avg.group.names
    # can we combine all those group types into a coherent graph? Need to send that data backwards up the chain as a data frame 
    # We know the returned results will have "avg" columns at the end of each but these are already converted to fold-change so it'll take some jury-rigging print instead
    # The data frame will look like $start.percent, $strain.name, $start.pop [all 1], $MX.RX.E0X.avg, ... and so on
    
    # Need to make a new group.data df to match up with the new column information - starts at 3, 3+(num.exp+1), 3+2*(num.exp+1)...
    new.col.list = seq(3, ncol(group.avg.df), num.exp+1)
    # rebuild the group information for a last set of line graphs
    group.avg.info = data.frame(rep.num =as.character(seq(1,length(new.col.list))), rep.name = c(paste(group.names, "avg", sep="."), "all.avg"), group = c(paste(as.character(group.names), "avg", sep="."), "all.avg"), stringsAsFactors=FALSE)
    
    output.directory = paste0(c("./", output.directory, "/rep.combined.avg/"), collapse="")
    dir.create(output.directory, recursive=TRUE)
    
    group.graph.data.list = generate.line.graph.data(group.avg.df, rownames(group.avg.df), exp.names, new.col.list, group.avg.info)   
    print.graph.data(rownames(group.avg.df), group.graph.data.list, group.avg.info, paste0(c(output.directory, "group.avg.GC"), collapse=""), length(unlist(unique(group.data$group))))
    
    group.graph.data.list = generate.line.graph.data(group.norm.avg.df, rownames(group.norm.avg.df), exp.names, new.col.list, group.avg.info)   
    print.graph.data(rownames(group.norm.avg.df), group.graph.data.list, group.avg.info, paste0(c(output.directory, "group.norm.avg.GC"), collapse=""), length(unlist(unique(group.data$group))))
    
  } 
  
  #make.rep.plot.data(data.df.rep, paste0(c(output.directory, "by.rep"), collapse=""), norm.strain, num.rep, num.exp, exp.names, group.file)
}


# Catch-all function to make the plots for information sheets for a given data file. 
# It will do all the work except you must format the file correctly and give the specific grouping columns.
# Other than that, it will create boxplots, correlation heatmaps, strain representation heatmaps and normalized versions of some of these
make.rep.plot.data = function(data.df, output.directory, norm.strain, num.rep, num.exp, exp.names, group.data) {
  
  exp.col = c(seq(2,ncol(data.df), by=floor(ncol(data.df)/num.rep)))
  
  #make.plot.data("140531.TF.M3.analysis/140531.M3.all.exp.txt", "140531.TF.M3.analysis/140606.test", "VC20019", c(1,2,3,13,29,45,55), c(1,2,3,4,6,8,10))
  
  # When you read the table in, if you don't specify a colClasses type, it looks like it will guess at the format. 
  # Fortunately it looks like it takes our data and makes it numeric. The row.names call makes the first column into our row.names and removes it from the table
  # data.df = read.table(data.file, header=TRUE, sep="\t", row.names=1)
  data.norm.df = as.data.frame(normalize.strains(data.df, norm.strain))
  
  #run.command = paste0(c("make.rep.plot.data(\"", data.file, "\", \"", output.directory, "\", \"", norm.strain, "\", c(", paste0(c(exp.col), collapse=","), "), c(", paste0(c(exp.names), collapse=","), "))"), collapse="")
  
  output.directory = paste0(c("./", output.directory, "/"), collapse="")
  # Create the output directory
  dir.create(output.directory, recursive=TRUE)
  
  #   output.info <<- paste0(c(output.directory, "run.info.txt"), collapse = "")
  #   
  #   lapply(c(run.command), write, output.info, append=FALSE, ncolumns=1000)
  
  lapply(c(paste0("Starting make.rep.plot.data at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)  
  
  # Steps we want to take
  # 1. generate correlations: both pearson and spearman
  # generate.cor.heatmap = function(pool.data.df, cor.method, output.file)
  
  heatmap.file = paste0(c(output.directory, "cor.heatmap.pearson"))
  data.cor.heat.pearson = generate.cor.heatmap(data.df, "pearson", heatmap.file)
  
  heatmap.file = paste0(c(output.directory, "cor.heatmap.spearman"))
  data.cor.heat.spearman = generate.cor.heatmap(data.df, "spearman", heatmap.file)
  
  # 2. generate heatplots of strain representation: both by relative abundance and normalized - this can be done in a single call
  strain.heatmap.file = paste0(c(output.directory, "strain.heatmap"))
  data.strain.heatmap = generate.strain.heatmap(data.df, norm.strain, strain.heatmap.file)
    
  # 5. Produce a series of line graphs representing the growth curves
  # do it for both abundance and normalized data
  # Currently the system will have to regenerate the growth curve data each time it's called (% vs. normalized)
  
  lapply(c("Generating growth curves and graphs for abundance data"), write, output.info, append=TRUE, ncolumns=1000)
  
  growth.curve.dir = paste0(c(output.directory, "growth.curve.data/"), collapse="")
  dir.create(growth.curve.dir)
  growth.curve.file = paste0(c(growth.curve.dir, "GC.abundance"), collapse="")
  data.growth.curve = growth.curve(data.df, 1, num.rep, num.exp, exp.names, growth.curve.file, group.data)
  #growth.curve = function(series.data, start.col, num.rep, num.exp, exp.names, output.file, group.file=NULL)
  #fold.change.analysis = function(strain.list, strain.graph.data.list, num.rep, num.exp, exp.names, output.file, group.data, norm.strain)
  
  # 6. Produce a series of heatmaps showing the fold-change and rate of fold-change both by strain and by rep
  # This will use the graph data produced in the last step
  
  strain.graph.data.list = data.growth.curve[[3]]
  
  fc.file.dir = paste0(c(output.directory, "heatmap.fold.change.data/"), collapse="")
  dir.create(fc.file.dir)
  fold.change.analysis(rownames(data.growth.curve[[1]]), strain.graph.data.list, num.rep, num.exp, exp.names, fc.file.dir, group.data, norm.strain, log2(data.df[,1]))
  
  ## ------------- 150323: NOTE on normalizing when norm.strain abundance = 0 ------------- ##
  # If the abundance of a normalizing strain in a particular column (rep) is 0, then you'll be generating a series of "Inf" values
  # This is especially troublesome when trying to correlate a matrix with infinite values
  
  lapply(c("Generating growth curves and graphs for normalized data"), write, output.info, append=TRUE, ncolumns=1000)
  
  growth.curve.file = paste0(c(growth.curve.dir, "GC.norm"), collapse="")
  data.norm.growth.curve = growth.curve(data.norm.df, 1, num.rep, num.exp, exp.names, growth.curve.file, group.data)
  #data.growth.curve = growth.curve(data.norm.df, 1, exp.col, exp.names, growth.curve.file, group.file)

  
  return (list(data.growth.curve, data.norm.growth.curve))
  
}

# Catch-all function to make the plots for information sheets for a given data file. 
# It will do all the work except you must format the file correctly and give the specific grouping columns.
# Other than that, it will creat boxplots, correlation heatmaps, strain representation heatmaps and normalized versions of some of these
make.exp.plot.data = function(data.df, output.directory, norm.strain, exp.col, exp.names, group.data, num.rep) {
  
  #make.plot.data("140531.TF.M3.analysis/140531.M3.all.exp.txt", "140531.TF.M3.analysis/140606.test", "VC20019", c(1,2,3,13,29,45,55), c(1,2,3,4,6,8,10))
  
  # When you read the table in, if you don't specify a colClasses type, it looks like it will guess at the format. 
  # Fortunately it looks like it takes our data and makes it numeric. The row.names call makes the first column into our row.names and removes it from the table
  #data.df = read.table(data.file, header=TRUE, sep="\t", row.names=1)
  data.norm.df = as.data.frame(normalize.strains(data.df, norm.strain))
  
  #run.command = paste0(c("make.exp.plot.data(\"", data.file, "\", \"", output.directory, "\", \"", norm.strain, "\", c(", paste0(c(exp.col), collapse=","), "), c(", paste0(c(exp.names), collapse=","), "))"), collapse="")
  
  output.directory = paste0(c(output.directory, "/"), collapse="")
  # Create the output directory
  dir.create(output.directory, recursive=TRUE)
  
  #output.info <<- paste0(c(output.directory, "run.info.txt"), collapse = "")
  
  lapply(c(paste0("Starting make.gen.plot.data at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  # Steps we want to take
  # 1. generate correlations: both pearson and spearman
  # generate.cor.heatmap = function(pool.data.df, cor.method, output.file)
  
  lapply(c(paste0("Generating correlation heat maps at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  heatmap.file = paste0(c(output.directory, "cor.heatmap.pearson"))
  data.cor.heat.pearson = generate.cor.heatmap(data.df, "pearson", heatmap.file)
  write.table(data.cor.heat.pearson[[2]], file=paste(c(heatmap.file, ".data.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  heatmap.file = paste0(c(output.directory, "cor.heatmap.spearman"))
  data.cor.heat.spearman = generate.cor.heatmap(data.df, "spearman", heatmap.file)
  write.table(data.cor.heat.spearman[[2]], file=paste(c(heatmap.file, ".data.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  # calculate the average correlation of the different generations
  num.gen = length(exp.names)-1
  data.avg.matrix = matrix(nrow=num.gen, ncol=num.gen)
  
  for (i in 1:num.gen) { 
    for (j in 1:num.gen) { 
      data.avg.matrix[i,j] = mean(as.numeric(unlist(data.cor.heat.spearman[[2]][(((i-1)*num.rep)+2):(i*num.rep), (((j-1)*num.rep)+2):(j*num.rep)])))
    }
  }
  
  heatmap.file = paste0(c(output.directory, "cor.heatmap.spearman.avg"))
  data.cor.heat.spearman.avg = generate.generic.heatmap(data.avg.matrix, heatmap.file)
  write.table(data.cor.heat.spearman.avg[[2]], file=paste(c(heatmap.file, ".data.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  
  # 2. generate heatplots of strain representation: both by relative abundance and normalized - this can be done in a single call
  lapply(c(paste0("Generating abundance heatmaps at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  strain.heatmap.file = paste0(c(output.directory, "strain.heatmap"))
  data.strain.heatmap = generate.strain.heatmap(data.df, norm.strain, strain.heatmap.file)
  
  # 3. generate boxplot charts by combining columns to represent specific expansions (or groups): both by relative abundance and normalized
  lapply(c(paste0("Generating boxplots at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  expansion.boxplot.dir = paste0(c(output.directory, "expansion.boxplot.data/"), collapse="")
  dir.create(expansion.boxplot.dir)
  
  expansion.outlier.df = data.frame()
  expansion.outlier.norm.df = data.frame()
  
  for (i in 1:(length(exp.col)-1)) {
    
    start.col = exp.col[i]
    end.col = exp.col[i+1] - 1
    
    # Note that you can't make a boxplot series on data sets of size 1. So we can only make them for expansions with at least 2 replicates
    if ((end.col-start.col) > 0) {
      expansion.boxplot.file = paste0(c(expansion.boxplot.dir, "generation.", exp.names[i]), collapse="")
      boxplot.info = make.boxplot.data(data.df, start.col, end.col, expansion.boxplot.file)
      
      # the outliers data won't have the current expansion name attached so we make one to go with it
      curr.outliers = boxplot.info[[3]]
      curr.outliers = cbind(rep(exp.names[i], nrow(curr.outliers)), curr.outliers)
      expansion.outlier.df = rbind(expansion.outlier.df, curr.outliers)
      
      # ----- Now do it with the normalized data ----- #
      
      expansion.boxplot.file = paste0(c(expansion.boxplot.dir, "generation.", exp.names[i], ".norm"), collapse="")
      boxplot.info = make.boxplot.data(data.norm.df, start.col, end.col, expansion.boxplot.file)  
      
      # the outliers data won't have the current expansion name attached so we make one to go with it
      curr.outliers = boxplot.info[[3]]
      curr.outliers = cbind(rep(exp.names[i], nrow(curr.outliers)), curr.outliers)
      expansion.outlier.norm.df = rbind(expansion.outlier.norm.df, curr.outliers)
      
    }
  }
  
  # Now do the last columns as long as there are more than 1
  start.col = tail(exp.col,1)
  end.col = ncol(data.df)
  
  if (end.col-start.col > 0) {
    
    expansion.boxplot.file = paste0(c(expansion.boxplot.dir, "generation.", tail(exp.names,1)), collapse="")
    boxplot.info = make.boxplot.data(data.df, tail(exp.col,1), ncol(data.df), expansion.boxplot.file)
    
    curr.outliers = boxplot.info[[3]]
    curr.outliers = cbind(rep(exp.names[i], nrow(curr.outliers)), curr.outliers)
    expansion.outlier.df = rbind(expansion.outlier.df, curr.outliers)
    
    # ----- Now do it with the normalized data ----- #
    
    expansion.boxplot.file = paste0(c(expansion.boxplot.dir, "generation.", tail(exp.names,1), ".norm"), collapse="")
    boxplot.info = make.boxplot.data(data.norm.df, tail(exp.col,1), ncol(data.norm.df), expansion.boxplot.file)
    
    curr.outliers = boxplot.info[[3]]
    curr.outliers = cbind(rep(exp.names[i], nrow(curr.outliers)), curr.outliers)
    expansion.outlier.norm.df = rbind(expansion.outlier.norm.df, curr.outliers)
    
  }
  
  # Now write it all to files
  write.table(expansion.outlier.df, file=paste(c(expansion.boxplot.dir, "outlier.data.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  write.table(expansion.outlier.norm.df, file=paste(c(expansion.boxplot.dir, "outlier.data.norm.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  # 4. generate growth curve as boxplot data by combining appropriate columns: both by relative abundance and normalize
  
  #     gc.dir = paste0(c(output.directory, "boxplot.growth.curve.abundance/"), collapse="")
  #     boxpot.info = make.boxplot.growth.curve(data.df, exp.col, exp.names, gc.dir, norm.strain="")
  #    
  #     gc.dir = paste0(c(output.directory, "boxplot.growth.curve.norm/"), collapse="")
  #     boxplot.info = make.boxplot.growth.curve(data.norm.df, exp.col, exp.names, gc.dir, norm.strain=norm.strain)
  
  # 4. generate heatplots of strain representation: both by relative abundance and normalized for every expansion
  
  num.reps = floor(ncol(data.df)/(length(exp.names)-1))
  exp.data.df.list = list()
  
  for (i in 1:(length(exp.names)-1)) {    
    
    lapply(c(paste0("picking dataset from ", 2 + num.reps*(i-1), " to ", num.reps*i)), write, output.info, append=TRUE, ncolumns=1000)
    
    strain.heatmap.file = paste0(c(output.directory, "strain.heatmap.generation.", exp.names[i+1]))
    exp.data.df = data.df[, c(1, (2 + num.reps*(i-1)):(1 + num.reps*i))]
    
    temp.data.df = data.df[, c((2 + num.reps*(i-1)):(1 + num.reps*i))]
    colnames(temp.data.df) = group.data$rep.name
    
    exp.data.df.list[[i]] = temp.data.df
    
    #data.strain.heatmap = generate.strain.heatmap(exp.data.df, norm.strain, strain.heatmap.file)
    data.strain.heatmap = generate.strain.heatmap(temp.data.df, norm.strain, strain.heatmap.file)
    
    # Generate individual correlation heatmaps
    heatmap.file = paste0(c(output.directory, "cor.heatmap.pearson.gen.", exp.names[i+1]))
    #data.cor.heat.pearson = generate.cor.heatmap(exp.data.df, "pearson", heatmap.file)
    data.cor.heat.pearson = generate.cor.heatmap(temp.data.df, "pearson", heatmap.file)
    
    heatmap.file = paste0(c(output.directory, "cor.heatmap.spearman.gen.", exp.names[i+1]))
    #data.cor.heat.spearman = generate.cor.heatmap(exp.data.df, "spearman", heatmap.file)
    data.cor.heat.spearman = generate.cor.heatmap(temp.data.df, "spearman", heatmap.file)
    
    
    
  }
  
  # 5. generate heatplots of each strain across the different expansions that they have been given (ie plot cols of reps, rows of expansions)
  # Added 150510
  
  strain.heatmap.dir = paste0(c(output.directory, "strain.gen.heatmaps.data/"), collapse="")
  dir.create(strain.heatmap.dir)
  
  strain.list = row.names(data.df)
  
  for (i in 1:length(strain.list)) {
    
    curr.strain = strain.list[i]
    curr.set.df = exp.data.df.list[[1]][0,]

    for (j in 1:length(exp.data.df.list)) {
      
      # Grab the same row (strain) and put it with each of the other expansion versions
      curr.set.df = rbind(curr.set.df, exp.data.df.list[[j]][i,])
      rownames(curr.set.df)[j] = paste0(c(curr.strain, ".gen.", exp.names[j+1]), collapse="")
    }
    
    strain.heatmap.file = paste0(c(strain.heatmap.dir, curr.strain, ".strain.heatmap.by.generation"))
    generate.strain.fold.change.heatmap(curr.set.df, strain.heatmap.file, "series", "generation", "abundance\n") 
  } 
  
}

# Given data that is likely in a data frame
# Convert it to a matrix, then normalizes it upon a given strain within the group
normalize.strains = function(data, norm.strain) {
  
  norm.strain.row = which(rownames(data)==norm.strain)
  
  # 151113: Add in a little cheat in case we have a 0 value for norm.strain. Otherwise we start getting weird and incorrect numbers
  #data.matrix[norm.strain.row[which(norm.strain.row == 0)] = 0.00001
  
  data[norm.strain.row, which(data[norm.strain.row,] == 0)] = 0.00001
  
  data.matrix = as.matrix(data)
    
  norm.data = t(t(data.matrix)/data.matrix[norm.strain.row,])

  return(norm.data)
}

# This takes the data and rebuilds it to be ready for a boxplot of a particular expansion
# relies on the reshape2 library
# Assumes the row.names are the strain names
# Does not work correctly on a matrix. data must be a data frame object
# This is really used to analyse pools of populations to give a quick picture of what overall representation is
# Altered 140610 to generate 2 kinds of plots: With and without outliers labeled as the labeled version can be very messy.
# We should also collect the outlier data into a file: Each row
make.boxplot.data = function(data, start.col, end.col, output.file) {
  
  #data = as.matrix(data)
  outlier.output.file = paste0(c(output.file, ".label.png"), collapse="")
  output.file = paste0(c(output.file, ".png"), collapse="")
  graph.title = tail(unlist(strsplit(output.file, "\\/")), 1)
  
  strain = rownames(data)
  data = cbind(strain, data[,start.col:end.col])
  new.data = melt(data, id.vars=c("strain"), var="series")
  outliers.df = data.frame()
  
  # For some reason boxplot doesn't set margins the same way so you need this command
  # see: https://stat.ethz.ch/pipermail/r-help/2009-April/196805.html
  par(mar=(c(6,5,2,1)))

  #----- First make the plot with the outliers unlabeled -----#
  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=14, 
      height=8.5, 
      pointsize=font.size, 
      res=300   )
  
  # The use of las=2 sets the x-labels to 90 degrees
  #boxplot.data = boxplot(value ~ reorder(strain, value), new.data, ylab="% Representation", las=2)
  # Remove outliers
  # boxplot.data = boxplot(value ~ reorder(strain, value), new.data, ylab="% Representation", las=2, outline=FALSE)
  # unordered with outliers
  boxplot.data = boxplot(value ~ strain, new.data, ylab="% Representation", las=2, outline=TRUE, main=graph.title, range=IQ.range)
  
  dev.off()
  
  #----- Now make the plot with the outliers labeled -----#
  
  png(filename=outlier.output.file, 
      type="cairo",
      units="in", 
      width=14, 
      height=8.5, 
      pointsize=font.size, 
      res=300   )
  
  # The use of las=2 sets the x-labels to 90 degrees
  #boxplot.data = boxplot(value ~ reorder(strain, value), new.data, ylab="% Representation", las=2)
  # Remove outliers
  # boxplot.data = boxplot(value ~ reorder(strain, value), new.data, ylab="% Representation", las=2, outline=FALSE)
  # unordered with outliers
  boxplot.data = boxplot(value ~ strain, new.data, ylab="% Representation", las=2, outline=TRUE, main=graph.title, range=IQ.range)
  
  # Can we make the boxplot interactive to expand only if there are outliers on the right side?
  if (length(boxplot.data$group[boxplot.data$group %in% length(strain)]) > 0) {
    
    # Then we need to replot the boxplot with a longer xlim to accomodate the text we'll be adding
    boxplot.data = boxplot(value ~ strain, new.data, ylab="% Representation", las=2, outline=TRUE, main=graph.title, range=IQ.range,
                           xlim=c(0.5, length(strain) + 3))
  } 
  
  if(length(boxplot.data$out > 0)) {
  
    for (i in 1:length(boxplot.data$out)) {
      
      outlier.data = new.data[which(new.data$strain == boxplot.data$names[boxplot.data$group[i]]),]
      curr.outlier = outlier.data[which(outlier.data$value == boxplot.data$out[i]),]
      outliers.df = rbind(outliers.df, curr.outlier)
      
      #text(boxplot.data$group[i], boxplot.data$out[i], outlier.data$series[which(outlier.data$value == boxplot.data$out[i])], pos=4)
      text(boxplot.data$group[i], boxplot.data$out[i], curr.outlier$series, pos=4)
    }
  
  }
  # unordered, no outliers
  #boxplot.data = boxplot(value ~ strain), new.data, ylab="% Representation", las=2, outline=FALSE)
  
  dev.off()
  
  # generate outliers as these could indicate something interesting like an interaction
  #outliers = new.data[new.data$value %in% boxplot.data$out,]
  
  return (list(new.data, boxplot.data, outliers.df))
}


# The purpose of this function is to generate "growth" curves given a series of data
# This essentially helps us to determine of a strain is growing or not
# A boxplot is generated from the data and saved to the output.directory based on the strain name
# the variable norm.strain is a strain name OR it is blank ("")
# Can't seem to melt the data properly so you'll have to make it by scratch.
make.boxplot.growth.curve = function(series.data, series.rows, group.list, output.directory, norm.strain) {
  
  
  #output.directory = paste(c("./", output.directory, "/"), collapse = "")
  # Create the output directory
  dir.create(output.directory)
  
  raw.data.list = list()
  boxplot.results = list()
  
  series.info = numeric()
  
  # Generate the expansion names
  #   for (j in 1:(length(series.rows)-1)) {
  #     
  #     start.row = series.rows[j]
  #     end.row = series.rows[j+1]
  #     
  #     series.info = c(series.info, rep(group.list[j], end.row-start.row))
  #   }
  # series.info = c(series.info, rep(tail(group.list,1), ncol(series.data)-tail(series.rows,1)+1))
  
  for (j in 1:(length(series.rows)-1)) {
    
    start.row = series.rows[j]
    end.row = series.rows[j+1]
    
    series.info = c(series.info, rep(j, end.row-start.row))
  }
  
  series.info = c(series.info, rep(length(series.rows), ncol(series.data)-tail(series.rows,1)+1))
  series.colnames = colnames(series.data)
  
  outliers.df = data.frame()

  
  # interate through each strain to make the boxplot data
  for (i in 1:nrow(series.data)) {
    
    curr.strain.outliers.df = data.frame()
    
    # get the strain name and then the data
    curr.strain = rownames(series.data[i,])
    
    lapply(c("Making boxplot growth curve for strain: ", curr.strain), write, output.info, append=TRUE, ncolumns=1000)
    
    # generate a filename
    curr.file.name = paste(c(output.directory, curr.strain, ".boxplot.growth.curve.png"), collapse = "")
    curr.file.name.clean = paste(c(output.directory, curr.strain, ".boxplot.growth.curve.clean.png"), collapse = "")
    
    # This will make a column matrix anyways
    curr.strain.data = as.matrix(as.numeric(series.data[i,]))
    
    curr.strain.data = cbind(series.info, curr.strain.data)
    
    colnames(curr.strain.data) = c("expansion", "value")
    
    #curr.strain.data.mirror = cbind(curr.strain.data, series.colnames)
    #colnames(curr.strain.data.mirror) = c("expansion", "value", "series")
    
    curr.strain.data = as.data.frame(curr.strain.data)
    #curr.strain.data.mirror = as.data.frame(curr.strain.data.mirror)
    
    # In the case of normalized strains, it can make the graphs hard to read so we'll make two versions: 
    # 1. with outliers removed
    # 2. with outliers (labeled)
    
    if (norm.strain !="") {
      
      # ----- 1. outliers removed ----- #
      
      par(mar=(c(6,5,2,4)))
            
      png(filename=curr.file.name.clean, 
          type="cairo",
          units="in", 
          width=11, 
          height=8.5, 
          pointsize=font.size, 
          res=300   )
      
      boxplot.data = boxplot(value ~ expansion, curr.strain.data, ylab=paste0(c("Abundance vs. ", norm.strain), collapse=""), xlab="Expansion", outline=FALSE, main = paste(c(curr.strain), " boxplot growth curve"), collapse="",
                             names=group.list, range=IQ.range)
      
      dev.off()
      
      # ----- 2. outliers included and labeled (more complicated to label) ----- #
      
      par(mar=(c(6,5,2,4)))
      
      png(filename=curr.file.name, 
          type="cairo",
          units="in", 
          width=11, 
          height=8.5, 
          pointsize=font.size, 
          res=300   )
      
      boxplot.data = boxplot(value ~ expansion, curr.strain.data, ylab=paste0(c("Abundance vs. ", norm.strain), collapse=""), xlab="Expansion", outline=TRUE, main = paste(c(curr.strain), " boxplot growth curve"), collapse="",
                             names=group.list, range=IQ.range)
            
      # Plot the outliers (hopefully)
      if (length(boxplot.data$out > 0)) {
        
        # Are there any outliers in the last group?
        if (length(boxplot.data$out[boxplot.data$group %in% length(series.rows)])) {
          
          # Redo the boxplot with the longer xlim value
          boxplot.data = boxplot(value ~ expansion, curr.strain.data, ylab=paste0(c("Abundance vs. ", norm.strain), collapse=""), xlab="Expansion", outline=TRUE, main = paste(c(curr.strain), " boxplot growth curve"), collapse="",
                                 xlim=c(0.5, length(group.list)+ 1.5), names=group.list, range=IQ.range)
          
        }
        
        
        for (i in 1:length(boxplot.data$out)) {
          
          # Here now our $expansion should values (which are not our label names) will be the same as those of $group
          outlier.data = curr.strain.data[which(curr.strain.data$expansion == boxplot.data$group[i]),]
          outlier.series = series.colnames[which(curr.strain.data$expansion == boxplot.data$group[i])]
          
          curr.outlier = outlier.series[which(outlier.data$value == boxplot.data$out[i])]
          curr.outlier = cbind(curr.outlier, outlier.data[which(outlier.data$value == boxplot.data$out[i]),])
          
          curr.strain.outliers.df = rbind(curr.strain.outliers.df, curr.outlier)
          
          #lapply(c("outlier data :", outlier.data), write, output.info, append=TRUE, ncolumns=1000)
          
          lapply(c("outlier :", outlier.series[which(outlier.data$value == boxplot.data$out[i])]), write, output.info, append=TRUE, ncolumns=1000)
          
          text(boxplot.data$group[i], boxplot.data$out[i], outlier.series[which(outlier.data$value == boxplot.data$out[i])], pos=4)
        }
      }
      
      dev.off()
      
    } else {
      
      par(mar=(c(6,5,2,4)))
      
      png(filename=curr.file.name, 
          type="cairo",
          units="in", 
          width=11, 
          height=8.5, 
          pointsize=font.size, 
          res=300   )
      
      curr.strain.data$value = curr.strain.data$value*100
      # outline = TRUE shows us outliers 
      boxplot.data = boxplot(value ~ expansion, curr.strain.data, ylab="% Representation", xlab="Expansion", outline=TRUE, main = paste(c(curr.strain), " boxplot growth curve"), collapse="",
                             names=group.list, range=IQ.range)  
      
      # Plot the outliers (hopefully)
      if (length(boxplot.data$out > 0)) {
        
        # Are there any outliers in the last group?
        if (length(boxplot.data$out[boxplot.data$group %in% length(series.rows)])) {
          
          # Redo the boxplot with the longer xlim value
          boxplot.data = boxplot(value ~ expansion, curr.strain.data, ylab="% Representation", xlab="Expansion", outline=TRUE, main = paste(c(curr.strain), " boxplot growth curve"), collapse="",
                                 xlim=c(0.5, length(group.list) + 1.5), names=group.list, range=IQ.range)  
          
        }
        
        for (i in 1:length(boxplot.data$out)) {
          
          # Here now our $expansion should values (which are not our label names) will be the same as those of $group
          outlier.data = curr.strain.data[which(curr.strain.data$expansion == boxplot.data$group[i]),]
          outlier.series = series.colnames[which(curr.strain.data$expansion == boxplot.data$group[i])]
          
          curr.outlier = outlier.series[which(outlier.data$value == boxplot.data$out[i])]
          curr.outlier = cbind(curr.outlier, outlier.data[which(outlier.data$value == boxplot.data$out[i]),])
          
          curr.strain.outliers.df = rbind(curr.strain.outliers.df, curr.outlier)
          
          #lapply(c("outlier data :", outlier.data), write, output.info, append=TRUE, ncolumns=1000)
          
          lapply(c("outlier :", outlier.series[which(outlier.data$value == boxplot.data$out[i])]), write, output.info, append=TRUE, ncolumns=1000)
          
          text(boxplot.data$group[i], boxplot.data$out[i], outlier.series[which(outlier.data$value == boxplot.data$out[i])], pos=4)
        }
      }
      
      dev.off()
    }
        
    curr.strain.outliers.df = cbind(rep(curr.strain, nrow(curr.strain.outliers.df)), curr.strain.outliers.df)
    outliers.df = rbind(outliers.df, curr.strain.outliers.df)
    
    raw.data.list[[i]] = curr.strain.data
    boxplot.results[[i]] = boxplot.data
    
  }
  
  # You can't rename the columns if the data frame is empty
  if (nrow(outliers.df) > 0) {
    
    colnames(outliers.df) = c("strain", "series", "group", "value")      
    write.table(outliers.df, file=paste(c(output.directory, "outlier.data.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  }
    
  return (list(raw.data.list, boxplot.results, outliers.df))
  
}

### ----------------------- Growth rate analysis of data -----------------------###

# Goal: Essentially want to generate data based on
# 1. normalizing against a strain within the group (usually VC20019)
# 2. normalizing against a start point for that strain

# 150102: What if we don't normalize against the start point for a strain? What does it's growth pattern look like? Can we compare it to VC20019?
# If we do not use "fold-change" as a scale, the shape of the graphs look exactly the same, just the y-axis scales are different.
# We should generate a final graph of 1) VC20019 fold-change 2) strain fold-change 3) strain fold-change normalized


# Why? So if we want to generate a growth graph with some knowledge of what is happening to a "normal" strain
# The first conversion gets us against the baseline strain while the second turns this into a ratio of where the strain started in the population
# We should try and generate this for each replicate as well as an amalgamated/averaged one 
# This will require 
# data.info: either a normalized or unnormalized version of the data. This way we can be flexible in how we generate it
# start.col: the starting expansion column to work off of as the original starting abundance
# group.list: the start point of each grouping. In most case this will be the start of a replicate with the expansions in order
# group.names: the list of names used to represent each group
# corr.range is a duple indicating the upper and lower range required for two things to be considered correlated (default is .75-1)
# You'll need to catch the cases where the start.exp value is 0
# group.file: $name, $group 
# data.growth.curve = growth.curve(data.df, 1, num.rep, num.exp, exp.names, growth.curve.file, group.data)

# 150327: We should try to generate correlation plots of each REP within a strain? Then we can see which ones are standing out!
# This could be a better way to view strain plots than line graphs?

growth.curve = function(series.data, start.col, num.rep, num.exp, exp.names, output.file, group.data) {
  
  # Generate the names of the strains in which the initial read for abundance is 0. These cannot generate proper growth curve data.
  no.gc.data = rownames(series.data[series.data[,start.col] == 0,])
  write.table(no.gc.data, file=paste0(c(output.file, "no.gc.analysis.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  # now only take the strains/series that do have appropriate values in the initial column
  series.data = series.data[series.data[,start.col] > 0,]
  
  # grab the start population
  start.pop = series.data[,start.col]
  # 150102 added to look at raw changes.
  #start.pop = start.pop/start.pop
  output.df = cbind(start.pop, row.names(series.data))
  output.list = list() # This will be similar to output.df except each expansion or "set" will be held as a list on its own
  colnames(output.df) = c("start.percent", "strain.name")
  rownames(output.df) = row.names(series.data)
  
  # Used to keep track of the growing/expanding group start points in output.df. This will be handy later!
  new.col.list = numeric()
  
  avg.df = 0
  
  # Plan, iterate through the data grabbing the number of appropriate columns, based on the number expansions and replicates
  # 
  
  for (i in 1:num.rep) {
        
    # Grab the start and end columns you'll be working with
 
    curr.start = 2+(i-1)*num.exp
    curr.end = curr.start + num.exp - 1
    
    # curr.df is a multi-column data frame of each expansion for a single rep
    curr.df = series.data[,curr.start:curr.end]
    curr.df = cbind(start.pop, curr.df)
    
    # Here you are normalizing against the start point to look at fold change.
    # 150518 altered formula so fold-change at time 0 = 0. Makes more sense right? Essentially just subtract 1 from calculation. Affects all downstream results.
    #curr.df = (curr.df / start.pop) - 1
    #curr.df = (curr.df / start.pop)
    
    # Log2 fold change can help with defining the biological significance of the data BUT when we hit 0, log2(0) == INF
    # We need to account for possible cases where we ARE hitting infinity and should we graph on log fold change or just fold change?
    # Our other option is to mirror the fold change by taking any results less than one and taking the inverse of it.
    curr.df = log2(curr.df / start.pop)
    
    output.df = cbind(output.df, curr.df)
    output.list[[i]] = curr.df
    
    # Stack the abundance data from each replicate for each expansion
    # avg. df ends up having as many columns as there are expansions (including expansion 1 timepoint)
    avg.df = avg.df + curr.df
    
    # Add the START column of the NEXT grouping
    #new.col.list = c(new.col.list, (ncol(output.df)+1))
    # You're always adding a certain number of columns to the new data frame 
    # technically it should be growing by num.exp+1, starting at number col 3 as there are two columns we've added prior to that
    
    new.col.list = c(new.col.list, (3+(i-1)*(num.exp+1)))
  }
  
  # copy the new.col.list for the normalized data version
  norm.col.list = new.col.list
  
  # We're about to add one more set of columns (average). Mark this in the new col.list as well
  new.col.list = c(new.col.list, (ncol(output.df)+1))
  avg.df = avg.df/num.rep
  
  # 161201: altered the average.df analysis to replace any zero entries with an extremely small number
  # this should negate the catching mechanism built later in the generate.line.graph.data algorithm that removes malformed data.
  # 180528: I think when we import the data now, this will never happen since we set the minimum value of a datapoint to 0.00054.
  avg.df[avg.df == 0] = 0.00001

  # Add the final average of that data
  output.list[[(length(output.list) + 1)]] = avg.df
  
  # update the column names of the average information
  colnames(avg.df) = unlist(lapply(exp.names, function(x) paste(c(x, ".avg"), collapse="")))
  
  # Make a copy of the output information without the average data. We will use this for normalization
  norm.output.df = output.df
  
  # Add the average data information for the un-normalized graph data
  output.df = cbind(output.df, avg.df)
 
  avg.df.print = cbind(rownames(avg.df), avg.df)
  colnames(avg.df.print)[1] = c("strain")
    
  write.table(avg.df.print, file=paste0(c(output.file, "all.avg.data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  # 180613 can we do an analysis on the average abundance information to identify different strain growth parameters? 
  # This would be helpful in identify poor and well-growing strains vs more normal (VC20019 strains)
  # I make these constructs every time I'm trying to do a norm analysis. Maybe I should just make a quick function or something?
  strain.avg.list = list()
  for (i in 1:4) {
    strain.avg.list[[i]] = data.frame(strain=as.character(avg.df.print[,1]), stringsAsFactors = FALSE)
    
  }
  strain.avg.list[[5]] = data.frame()
  
  # We want to skip column 1 (strain names) and column 2 (start.population average = 0.00001)
  for (i in 3:ncol(avg.df.print)) {
    
    curr.avg = avg.df.print[,i]
    curr.avg.name = paste0(c("gen", exp.names[i-1], "avg"), collapse=".")
    curr.output = paste0(c(output.file, "strain.avg.fc.analysis", curr.avg.name), collapse=".")

    strain.avg.list[[1]][curr.avg.name] = curr.avg
      
    strain.avg.list = generate.norm.analysis(curr.avg, "\n average fold-change per strain\n Log2", curr.avg.name, strain.avg.list, curr.output, percent.outliers = 0.4, group.data=NULL)
    
  }
  
  generate.norm.files(strain.avg.list, output.file, "strain.avg.fold-change", xlab="Generation")
  
  write.table(output.df, file=paste0(c(output.file, "all.data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  # Now that we've generated an output file, we should make try to analyse the data quickly by
  # 1: Generating heatmaps of each set/replicate to see which strains look most like each other
  # 2: Generating line graphs of each strain as a represenation of all those strains
  # 3: Generating amalgamated line graphs (avg) for strains that are highly correlated with each other?
  
  
  # 1: Pearson correlation heatmaps by strain for each of the replicates and graph those strains that look the same together
  
  #   lapply(c("growth.curve: Making growth curves and heatmaps for entire data set"), write, output.info, append=TRUE, ncolumns=1000)
  #   
  #   for (i in 1:length(output.list)) {
  #     
  #     curr.set = output.list[[i]]
  #     curr.rep.num = group.data$rep.num[i]
  #     
  #     write.table(curr.set, file=paste(c(output.file, "rep", curr.rep.num, "txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=TRUE, col.names=TRUE, quote=FALSE)
  #     
  #     curr.set = t(as.matrix(curr.set))
  #     
  #     #lapply(c("Making growth curve for set ", i), write, output.info, append=TRUE, ncolumns=1000)
  #     
  #     if (i == length(output.list)) {
  #       curr.file.name = paste(c(output.file, "strain.cor.rep.avg"), collapse=".")
  #       
  #       # Need to do something after this to analyse the average and find heatmaps
  #       generate.cor.heatmap(curr.set, c("pearson"), curr.file.name)
  #       
  #       curr.file.name = paste(c(output.file, "strain.cor.rep.spearman.avg"), collapse=".")
  #       generate.cor.heatmap(curr.set, c("spearman"), curr.file.name)
  #       
  #       # Parameters: Correlation score 0.75-1      
  #       
  #     } else {
  #       curr.file.name = paste(c(output.file, "strain.cor.rep", curr.rep.num), collapse=".")
  #       generate.cor.heatmap(curr.set, c("pearson"), curr.file.name)
  #       
  #       curr.file.name = paste(c(output.file, "strain.cor.rep.spearman", curr.rep.num), collapse=".")
  #       generate.cor.heatmap(curr.set, c("spearman"), curr.file.name)
  #     }    
  #   }
  
  # 2: Generate line graphs for each strain
  
  strain.list = rownames(output.df)
  #lapply(c("Making graph data. new.col.list is ", new.col.list), write, output.info, append=TRUE, ncolumns=1000)
  
  # For some reason, output.df is a matrix at this point so we cast it back to a data frame for future manipulations
  output.df = as.data.frame(output.df)

  # The list that is returned is a list of data.frames() where the colnames are of the dfs are (strain.entry, type, [exp.names])
  lapply(c("growth.curve: Starting to generate line graph data"), write, output.info, append=TRUE, ncolumns=1000)
  strain.graph.data.list = generate.line.graph.data(output.df, strain.list, exp.names, new.col.list, group.data)

  # We need to make a graph with the average for each strain and the average population growth?
  avg.df.graph = cbind(avg.df.print[,1], rep("strains", length(strain.list)), avg.df.print[,2:ncol(avg.df.print)])
  
  
  # This line is not adding the right type of information
  colnames(avg.df.graph) = c("strain.entry", "type", exp.names)
  #temp.df = list(c("Popn", "pop.mean", apply(avg.df.graph[,3:ncol(avg.df.graph)],2,mean)))
  levels(avg.df.graph$strain.entry) = c(as.character(avg.df.graph$strain.entry), "Popn")
  levels(avg.df.graph$type) = c(as.character(avg.df.graph$type), "pop.mean")
  avg.df.graph = rbind(avg.df.graph, c("Popn", "pop.mean", apply(avg.df.graph[,3:ncol(avg.df.graph)],2,mean)))
  #avg.df.graph = rbind(avg.df.graph, temp.df)

  write.table(avg.df.graph, file=paste0(c(output.file, "population.avg.line.graph.data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

  # Need to make a quick temporary variable to print up all the graphs including the average population information
  strain.graph.data.list.temp = strain.graph.data.list
  strain.graph.data.list.temp[[length(strain.list)+1]] = avg.df.graph
  strain.list.print = c(strain.list, "Popn.mean")
  
  # Print out the line graph data for use later. Maybe for making quick adjustments.
  for (i in 1:length(strain.graph.data.list.temp)) {    
    write.table(strain.graph.data.list.temp[[i]], file=paste0(c(output.file, strain.list.print[i], "line.graph.data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }

  lapply(c("growth.curve: Starting to print line graph data"), write, output.info, append=TRUE, ncolumns=1000)
  #print.graph.data(c(strain.list, "pop.mean"), strain.graph.data.list.temp, group.data, output.file, y.axis="fold-change (Log2)")
  print.graph.data(c(strain.list, "pop.mean"), strain.graph.data.list.temp, group.data, output.file, y.axis=expression('fold-change (log'[2]*')'))


  #3a Generate normalized by division of avg line graphs for each strain

  # 151112: Normalize the data a second time by dividing the average into each set
  
  div.output.df = norm.output.df
  
  for (i in 1:num.rep) {
    
    # Grab the start and end columns you'll be working with
    
    #curr.start = 2+(i-1)*num.exp
    #curr.end = curr.start + num.exp - 1
    curr.start = norm.col.list[i]
    curr.end = curr.start + num.exp
    
    div.output.df[,curr.start:curr.end] = div.output.df[,curr.start:curr.end]/avg.df
    
  }
  
  write.table(div.output.df, file=paste0(c(output.file, "all.data.norm.div.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  
  div.output.df = as.data.frame(div.output.df)
  
  # The list that is returned is a list of data.frames() where the colnames are of the dfs are (strain.entry, type, [exp.names])
  lapply(c("growth.curve: Starting to generate normalized line graph data"), write, output.info, append=TRUE, ncolumns=1000)
  strain.graph.data.div.list = generate.line.graph.data(div.output.df, strain.list, exp.names, norm.col.list, group.data)
  
  # Print out the line graph data for use later. Maybe for making quick adjustments.
  for (i in 1:length(strain.graph.data.div.list)) {    
    write.table(strain.graph.data.div.list[[i]], file=paste0(c(output.file, strain.list[i], "line.graph.data.norm.div.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
  
  lapply(c("growth.curve: Starting to print normalized line graph data"), write, output.info, append=TRUE, ncolumns=1000)
  div.output.file = paste0(c(output.file, "avg.norm.div"), collapse=".")
  print.graph.data(strain.list, strain.graph.data.div.list, group.data, div.output.file, y.axis="Fold-change (Percent of average)")
  
  #3b Generate normalized by subtraction of avg line graphs for each strain
  
  # 151112: Normalize the data a second time by subtracting the average from each set
  
  diff.output.df = norm.output.df
  
  for (i in 1:num.rep) {
    
    # Grab the start and end columns you'll be working with
    
    #curr.start = 2+(i-1)*num.exp
    #curr.end = curr.start + num.exp - 1
    curr.start = norm.col.list[i]
    curr.end = curr.start + num.exp
    
    diff.output.df[,curr.start:curr.end] = diff.output.df[,curr.start:curr.end]-avg.df
    
  }
  
  write.table(diff.output.df, file=paste0(c(output.file, "all.data.norm.diff.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  
  diff.output.df = as.data.frame(diff.output.df)
  
  # The list that is returned is a list of data.frames() where the colnames are of the dfs are (strain.entry, type, [exp.names])
  lapply(c("growth.curve: Starting to generate normalized differential line graph data"), write, output.info, append=TRUE, ncolumns=1000)
  strain.graph.data.diff.list = generate.line.graph.data(diff.output.df, strain.list, exp.names, norm.col.list, group.data)
  
  # Print out the line graph data for use later. Maybe for making quick adjustments.
  for (i in 1:length(strain.graph.data.diff.list)) {    
    write.table(strain.graph.data.diff.list[[i]], file=paste0(c(output.file, strain.list[i], "line.graph.data.norm.diff.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
  
  lapply(c("growth.curve: Starting to print normalized line graph data"), write, output.info, append=TRUE, ncolumns=1000)
  diff.output.file = paste0(c(output.file, "avg.norm.diff"), collapse=".")
  print.graph.data(strain.list, strain.graph.data.diff.list, group.data, diff.output.file, y.axis="Fold-change (Difference from average")
  
  
  
  
  # 3: Generate correlation heatmaps of individual strains to compare their growth curves against each other per RNAi
  # Added this function on 150327
  
  #   for (i in 1:length(strain.graph.data.list)) {
  #     
  #     curr.strain = strain.list[i]
  #     # Grab the graph data
  #     curr.graph.data = strain.graph.data.list[[i]]
  #     # Name the rows and then prepare to cut the non-numeric entries
  #     rownames(curr.graph.data) = curr.graph.data$strain.entry
  #     # Keep the numeric data
  #     curr.graph.data = curr.graph.data[,3:ncol(curr.graph.data)]
  #     # Send it to be correlated
  #     curr.graph.data = t(as.matrix(curr.graph.data))
  #     
  #     
  #     curr.file.name = paste(c(output.file, curr.strain, "cor.by.rep"), collapse=".")
  #     generate.cor.heatmap(curr.graph.data, c("pearson"), curr.file.name)
  #     
  #     curr.file.name = paste(c(output.file, curr.strain, "cor.by.rep.spearman"), collapse=".")
  #     generate.cor.heatmap(curr.graph.data, c("spearman"), curr.file.name)
  #     
  #   }
  
#   fc.heatmap.data.list = list()
#     
#   # 4: Generate heatmaps versions of the line graph data 
#   for (i in 1:length(strain.graph.data.list)) {    
#     
#     curr.strain = strain.list[[i]]
#     curr.graph.df = strain.graph.data.list[[i]]
#     
#     # Format the graph data by keeping only the values
#     rownames(curr.graph.df) = curr.graph.df[,1]
#     curr.graph.df = curr.graph.df[,3:ncol(curr.graph.df)]
#     curr.graph.matrix = as.matrix(curr.graph.df)
#     
#     colnames(curr.graph.matrix) = unlist(lapply(colnames(curr.graph.matrix), function(x) paste0(c(curr.strain, "gen", x), collapse=".")))
#     
#     # Data is formatted as rows of conditions and columns as generations
#     fc.heatmap.data.list[[i]] = curr.graph.matrix
#     
#     strain.heatmap.file = paste0(c(output.file, curr.strain, "strain.heatmap.by.fold-change"), collapse=".")
#     generate.strain.fold.change.heatmap(t(curr.graph.matrix), strain.heatmap.file, "series", "generation", "fold-change\n")
#     
#     # consider writing the heatmap data file to disk?
#     
#   }
# 
#   # 5: Generate heatmaps of fold-change data on a per-replicate basis. 
#   # We need to rifle through fc.heatmap.data.list and re-assemble from X heatmaps (num strains) to Y heatmaps (num replicates)
#   # One heatmap for each replicate ordered by strain (rows)
#   
#   fc.heatmap.data.by.rep.list = list()
#   
#   # First we have to generate place holders for each matrix in the list
#   curr.graph.matrix = fc.heatmap.data.list[[1]]
#   colnames(curr.graph.matrix) = unlist(lapply(exp.names, function(x) paste0(c("gen", x), collapse=".")))
#   #rep.names = rownames(curr.graph.matrix)
#   
#   
#   for (i in 1:nrow(curr.graph.matrix)) {
#     fc.heatmap.data.by.rep.list[[i]] = curr.graph.matrix[0,]
#   }
#   
#   for (i in 1:length(fc.heatmap.data.list)) {
#     
#     curr.heatmap = fc.heatmap.data.list[[i]]
#     colnames(curr.heatmap) = unlist(lapply(exp.names, function(x) paste0(c("gen", x), collapse=".")))
#     
#     # At the end of this you'll have distributed one row from each heatmap to it's appropriate replicated
#     for (j in 1:nrow(curr.heatmap)) {
#       
#       #curr.rep.heatmap = fc.heatmap.data.by.rep.list[[j]]
#       
#       #curr.rep.heatmap = rbind(curr.rep.heatmap, curr.heatmap[j,])
#       
#       fc.heatmap.data.by.rep.list[[j]] = rbind(fc.heatmap.data.by.rep.list[[j]], curr.heatmap[j,])
#       #fc.heatmap.data.by.rep.list[[j]] = curr.rep.heatmap
#       
#       if (i == length(fc.heatmap.data.list)) {
#         # You're on the last set so now you can name all the rows
#         rownames(fc.heatmap.data.by.rep.list[[j]]) = unlist(strain.list)
#         # Now you should try to plot all this data?
#         
#         heatmap.file = paste0(c(output.file, rownames(curr.heatmap)[j], "rep.heatmap.by.fold-change"), collapse=".")
#         generate.strain.fold.change.heatmap(t(fc.heatmap.data.by.rep.list[[j]]), heatmap.file, "generation", "strain", "fold-change\n")
#       }      
#     }
#   }

  # 6: Can we generate a fold-change comparison? What is the rate of change at each timepoint?
  # Does this normalize our data a little better?
  # Problems? If VC20019 is on a downward trajectory and the strain to compare is going up, then we'll get a large negative
  # If conversely sloping we'll get a very small negative value
  # Note: We can't do this with normalized data
  # Try to normalize data first by finding the VC20019 replicates and then applying to all the other strains?
  # fc.heatmap.data.list
  # strain.list
  
  

  # 4: Generate heatmap comparison using statmod library and compareGrowthCurves() function
  # Added this function on 150327
  
#   for (i in 1:length(strain.graph.data.list)) {
#     
#     curr.strain = strain.list[i]
#     # Grab the graph data
#     curr.graph.data = strain.graph.data.list[[i]]
#     # Name the rows and then prepare to cut the non-numeric entries
#     rownames(curr.graph.data) = curr.graph.data$strain.entry
#     # Keep the numeric data
#     curr.graph.data = curr.graph.data[,3:ncol(curr.graph.data)]
#     # Send it to be correlated
# #     curr.GC.comparison = compareGrowthCurves(rownames(curr.graph.data), curr.graph.data, nsim=1000)
# #     
# #     curr.file.name = paste(c(output.file, curr.strain, "compare.GC.txt"), collapse=".")
# #     write.table(curr.GC.comparison, file=curr.file.name, append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
# #     
#     
#     curr.file.name = paste(c(output.file, curr.strain, "cor.by.rep"), collapse=".")
#     generate.cor.heatmap(curr.graph.data, c("pearson"), curr.file.name)
#     
#     curr.file.name = paste(c(output.file, curr.strain, "cor.by.rep.spearman"), collapse=".")
#     generate.cor.heatmap(curr.graph.data, c("spearman"), curr.file.name)
#     
#   }

  return(list(output.df, output.list, strain.graph.data.list))
  
}

# helper function used to generate -at least- the table for a strain (or set of strains) that can be made into a line graph
# Goes through data.df and pulls out each strain (should be a single row) that is formatted column-wise to hold the replicate data
# col.names should be the final names of the data-points as columns of $replicate.data
# Then the columns are dissected into a single replicate per row $strain $type $replicate.data
# group.file will be formatted as $ rep.num, $rep.name, $group 
# Right now the line.graph.data is generated as either replicate data or as the avg data line. We want the ability to assign groups of replicates to
# a specific type, perhaps based on the growth media, conditions, or specific RNAi types.

# 151116 Update: Generate additional lines to highlight 2-standard deviations from the average.
# If possible, going to only highlight/name the lines that are outside of this area.

generate.line.graph.data = function(data.df, strain.list, col.names, set.col, group.data) {
  
  new.strain.data = list()
  # Going to add this extra end point in order to make life easier with the for loops 
  #set.col = c(set.col, ncol(data.df) + 1)
  
  for (i in 1:length(strain.list)) {
    
    # grab the current strain, then the data for that strain - it should only be a single row
    curr.strain = strain.list[i]
    curr.strain.data = data.df[rownames(data.df) %in% curr.strain,]
    
    curr.strain.reformat = data.frame()
    
    lapply(c(paste0("curr.strain is ", curr.strain)), write, output.info, append=TRUE, ncolumns=1000)
    #lapply(c(paste0("ncol curr.strain.data is ", ncol(curr.strain.data))), write, output.info, append=TRUE, ncolumns=1000)
    
    
    # Break down the data into each replicate
    for (j in 1:length(set.col)) {
      
      curr.set = data.frame()
      #lapply(c(paste0("curr col set is ", j)), write, output.info, append=TRUE, ncolumns=1000)
      
      # reformat the rows to match the col.name information      
      if (j == length(set.col)) {
        
        # The last set is always the average data from the group of data
        set.start = set.col[j]
        set.end = ncol(curr.strain.data)
        
        # grab the sub-columns you want - still a single row
        curr.set = curr.strain.data[, set.start:set.end]
        
        if (length(set.col) > nrow(group.data)) {
          
          rownames(curr.set) = paste0(c(curr.strain, "avg"), collapse=".")
          curr.set = cbind("avg", curr.set)
          
        } else {
          
          #rownames(curr.set) = paste0(c(curr.strain, group.data$rep.num[j]), collapse=".")
          rownames(curr.set) = group.data$rep.name[j]
          curr.set = cbind(group.data$group[j], curr.set)
        }
        
        
      } else { 
        
        # We are going through the row of data replicate by replicate
        set.start = set.col[j]
        set.end = set.col[j+1] - 1
        
        # grab the sub-columns you want
        curr.set = curr.strain.data[, set.start:set.end]
        
        # ----------------------- SET THE LINE GRAPH LABELS ------------------------------ #
        # Now we name the row based on information in the group.data$rep.name column
        #rownames(curr.set) = paste0(c(group.data$group[j], group.data$rep.num[j]), collapse=".")
        rownames(curr.set) = group.data$rep.name[j]
        
        # Now we need to give this to the correct group for later analysis

        curr.set = cbind(group.data$group[j], curr.set)
      }
      
      # 150323: We only want to keep the rows that have good data in them
      # 161201: This is causing an issue! What we should do instead is prevent any of the average line values from being 0.
      # An average line with a value of 0 means all values are 0, so we can substitute with something else and the division should be fine.
      if (!(any(curr.set == Inf) || (any(is.na(curr.set))))) {
        
        colnames(curr.set) = c("type", col.names)
        # bind the rows together      
        curr.strain.reformat = rbind(curr.strain.reformat, curr.set)     
        
      } else {
        lapply(c(paste0("removing rep ", j)), write, output.info, append=TRUE, ncolumns=1000)
      }
       
    }    
    
    # Move the rownames into a column you can access for making the graph data (This comes from group.data$rep.name)
    curr.strain.reformat = cbind(rownames(curr.strain.reformat), curr.strain.reformat)
        
    # This is the correct renaming of 
    colnames(curr.strain.reformat) = c("strain.entry", "type", col.names)
    
    
    # 151117 add standard deviation lines and curves
    # grab all the lines prior and generate a standard deviation for each column but don't use the average in the calculation
    start.col = 3
    end.col = ncol(curr.strain.reformat)
    num.sigma = 2
    
    # calculate the data
    stdev = apply(curr.strain.reformat[1:(nrow(curr.strain.reformat)-1),start.col:end.col], 2, sd)
    avg.data = curr.strain.reformat[nrow(curr.strain.reformat), start.col:end.col]
    upper.stdev = avg.data + num.sigma*stdev
    lower.stdev = avg.data - num.sigma*stdev
    
    # add it into the graph
    #curr.strain.reformat[nrow(curr.strain.reformat)+1,] = c("*", paste0(c(num.sigma, "sigma"), collapse="-"), upper.stdev)
    #curr.strain.reformat[nrow(curr.strain.reformat)+1,] = c("x", paste0(c(num.sigma, "sigma"), collapse="-"), lower.stdev)
    #curr.strain.reformat[nrow(curr.strain.reformat)+1,] = c("test", "2-sigma", upper.stdev)
    #curr.strain.reformat[nrow(curr.strain.reformat)+1,] = c("test2", "2-sigma", lower.stdev)
    
    #stdev.data = rbind(c("*", paste0(c(num.sigma, "sigma"), collapse="-"), as.numeric(upper.stdev)), c("x", paste0(c(num.sigma, "sigma"), collapse="-"), as.numeric(lower.stdev)))
    #colnames(stdev.data) = c("strain.entry", "type", col.names)
    #curr.strain.reformat = rbind(curr.strain.reformat, stdev.data)
    
    # 150323 Addition: What happens if data is sparse or for some reason, we get Inf and NA data? 
    # If we remove it at this point, the avg.values are already corrupted by the addition of NA values.
    # Try and catch it in earlier steps but clean up the data graphs here?
    
    new.strain.data[[i]] = curr.strain.reformat
        
  }
  
  names(new.strain.data) = strain.list
  
  return (new.strain.data)
}

print.graph.data = function(strain.list, strain.graph.data.list, group.data=NULL, output.file, num.groups=NULL, y.axis="Total abundance (fold-change)") {
  
  for (i in 1:length(strain.list)) {
    
    lapply(c(paste0("print.graph.data: curr.strain is ", strain.list[i])), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.strain.data = strain.graph.data.list[i]
    curr.strain.melt = melt(curr.strain.data, id.vars=c("strain.entry", "type"), variable.name="generation", value.name="fold.change")
    curr.strain.melt$generation = as.double(as.character(curr.strain.melt$generation))
    curr.strain.melt$fold.change = as.double(as.character(curr.strain.melt$fold.change))
    
    curr.strain.file = paste0(c(output.file, strain.list[i], "all.curves.png"), collapse=".")
    
    png(filename=curr.strain.file, 
        type="cairo",
        units="in", 
        width=12, 
        height=10, 
        pointsize=font.size, 
        res=300   )
    
    # Colorblind friendly palettes that don't necessarily work when projected
    #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #cbbPalette <- c("#000000", "#D55E00", "#0072B2", "#006600", "#FF33CC", "#9933FF")
    cbbPalette <- c("#000000", "#D55E00", "#0072B2", "#006600", "#FF33CC", "#9933FF", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#990000")
    
    # Needed to add this so the final combined averages can be consistent with the "average" still be solid red in colour
    # Also need to include an extra set of standard dev lines
    if (is.null(num.groups)) {      
      num.groups = length(unlist(unique(group.data$group)))
    }
  
    
    # 150113: Played with this code for 2+ days to try and get the legend to display both color and geom_point shape without any success.
    # The colour settins appear to override the line-point shapes with these large circles. 
    
    plot <- (ggplot(data=curr.strain.melt, aes(x=generation, y=fold.change, group=strain.entry, shape=as.factor(type), linetype=as.factor(type), colour=as.factor(type))) +
               # remove the grey background
               theme_classic() +
               #theme_bw() +
               # Set the "theme" of the plot including font size for labels and legends            
               theme(axis.text=element_text(size=font.size), axis.title=element_text(size=font.size), plot.title=element_text(size=font.size), 
                    legend.title=element_text(size=font.size), legend.text=element_text(size=font.size),
                    panel.grid.major.y = element_line(size=0.25, linetype = 'solid', colour = "grey"), panel.grid.minor = element_blank()) +
               # Add a line but not too thick
               geom_line(size=0.75) +
               # Add points based on the type or grouping
               #geom_point(aes(shape=type), size=3) + 
               geom_point(size=4) +        
               #scale_shape_discrete(name="Experimental\nCondition") +
               #scale_size_manual(values=c(1,10)) +
               # The linetype will be set by how many groups there are
               #scale_linetype_manual(name="Experimental\nCondition", values=c(2:(1+num.groups), 1), labels=sort(unlist(unique(curr.strain.melt$type)))) +
               scale_linetype_manual(values=c(3:(2+num.groups), 1, 2), guide=FALSE) +
               #scale_shape_manual(aes(colour=type, shape=type),name="Experimental\nCondition", values=c(16:(15+num.groups), 15), labels=sort(unlist(unique(curr.strain.melt$type)))) + 
               #scale_shape_manual(name="Experimental\nCondition", values=c(16:(15+num.groups), 15), labels=sort(unlist(unique(curr.strain.melt$type)))) + 
               #scale_shape_manual(values=c(17:(16+num.groups), 15, 16)) + 
               scale_shape_manual(values=c(0:(num.groups-1), 15, 16)) + 
               
               # Pick colours for the lines based on the adjusted color-blind friendly palette
               #scale_colour_manual(name="Experimental\nCondition", values=c(cbbPalette[1:num.groups], "#CC0033"), labels=sort(unlist(unique(curr.strain.melt$type))), guide=FALSE) +
               scale_colour_manual(values=c(cbbPalette[1:num.groups], "#CC0033", "#CC0033")) +
               labs(colour="Experimental\nGroup", shape="Experimental\nGroup", linetype="Experimental\nGroup") +
               #guides(colour="legend", size="none", shape="legend") +
               
               # ---- TITLES, LABELS and AXIS SETTINGS ---- #
               ggtitle(paste0(c("Growth curve for ", strain.list[i]), collapse="")) +
               
               # Use the direct.label library and the last.points function. Push the text just to the right of the end point with hjust but you need an expanded X-axis
               # Also when you use this line, you will end up making the legend incorrect - it will have the correct colour but will not put the correct geom_point shape.
               #geom_dl(aes(label=strain.entry, colour=type), method=list("last.points", hjust=-0.1)) + 
               #geom_dl(aes(label=strain.entry, colour=type), method=list("last.bumpup", cex=0.8, hjust=-0.1)) + 
               #scale_x_continuous(breaks=seq(min(curr.strain.melt$generation), max(curr.strain.melt$generation), 1), limits=c(min(curr.strain.melt$generation),max(curr.strain.melt$generation)+2)) +
               scale_x_continuous(breaks=seq(min(curr.strain.melt$generation), max(curr.strain.melt$generation), 1), limits=c(min(curr.strain.melt$generation),max(curr.strain.melt$generation))) +
               
               
               #scale_x_continuous(breaks=seq(min(curr.strain.melt$generation), max(curr.strain.melt$generation), 1), limits=c(min(curr.strain.melt$generation),max(curr.strain.melt$generation)+1.75)) +
               #scale_x_continuous(breaks=seq(min(curr.strain.melt$expansion), max(curr.strain.melt$expansion), 1), limits=c(1,max(curr.strain.melt$expansion))) +
               # Set the legend name
               #scale_shape_discrete(name="Experimental\nCondition") +
               #scale_color_discrete(name="Experimental\nCondition") +
               scale_y_continuous(name=y.axis, limits=c(min(0, min(curr.strain.melt$fold.change)),max(curr.strain.melt$fold.change))) 
               # We set the x-limit to +n to accomodate the direct labels and keep them on the graph
             
    )
    

    #print.graph.data(c("VC20019"), test.graph.data.list, group.data, "140113.graph.test", 3)
    #generate.plots("./141218.M10.E02.re-analysis.2/M12-E01.M11.E02-E08.data.by.exp.txt", "141218.M10.E02.re-analysis.2/150112.M11-E08.analysis.labels.v3", 24,4, "./141218.M10.E02.re-analysis.2/M10.M11.rep.info.2.txt") 
    print (plot)
    
    dev.off()
   
  }
  
}

# 150516: New function to do all the strain heatmap analysis that is analagous to the growth curves
# 151118: this function only gets used on the growth curve data sent back from the growth curve function earlier.
# The growth.curve generation now adds 2 more rows (columns on the heatmap) of data that represent 2-sigma.
# This data in particular is not very useful when viewing as a heatmap.
# 180426: I think the 2-sigma data was removed from the growth curve data. Looks like it's just the reps and an avg row now. 

fold.change.analysis = function(strain.list, strain.graph.data.list, num.rep, num.exp, exp.names, output.file, group.data, norm.strain, log.start.pop) {
  
  fc.heatmap.data.list = list()
  
  # 4: Generate heatmaps versions of the line graph data 
  for (i in 1:length(strain.graph.data.list)) {    
    
    curr.strain = strain.list[[i]]
    curr.graph.df = strain.graph.data.list[[i]]
    
    # Format the graph data by keeping only the values
    rownames(curr.graph.df) = curr.graph.df[,1]
    # Below line cuts off the avg and "2-sigma" rows. 180426, change the row subtraction to 1 since only average row is present.
    curr.graph.df = curr.graph.df[1:(nrow(curr.graph.df)-1),3:ncol(curr.graph.df)]
    curr.graph.matrix = data.matrix(curr.graph.df)
    
    colnames(curr.graph.matrix) = unlist(lapply(colnames(curr.graph.matrix), function(x) paste0(c(curr.strain, "gen", x), collapse=".")))
    
    # Data is formatted as rows of conditions and columns as generations
    fc.heatmap.data.list[[i]] = curr.graph.matrix
    
    strain.heatmap.file = paste0(c(output.file, curr.strain, ".strain.heatmap.by.fold-change"), collapse="")
    generate.strain.fold.change.heatmap(t(curr.graph.matrix), strain.heatmap.file, "series", "generation", "fold-change\n")
    
    curr.graph.matrix = cbind(rownames(curr.graph.matrix), curr.graph.matrix)
    colnames(curr.graph.matrix)[1] = "series"
    # This data file is nearly the same as the line graph data in terms of numbers but some columns are changed and the averaged out row at the bottom is removed.
    write.table(curr.graph.matrix, file=paste0(c(strain.heatmap.file, "data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

  }
  
  # 5: Generate heatmaps of fold-change data on a per-replicate basis. 
  # We need to rifle through fc.heatmap.data.list and re-assemble from X heatmaps (num strains) to Y heatmaps (num replicates)
  # One heatmap for each replicate ordered by strain (rows)
  
  fc.heatmap.data.by.rep.list = list()
  
  # First we have to generate place holders for each matrix in the list
  curr.graph.matrix = fc.heatmap.data.list[[1]]
  colnames(curr.graph.matrix) = unlist(lapply(exp.names, function(x) paste0(c("gen", x), collapse=".")))
  #rep.names = rownames(curr.graph.matrix)
  
  
  # Set the number column names for the list of matrices you want to generate n rows represent n strains?
  for (i in 1:nrow(curr.graph.matrix)) {
    fc.heatmap.data.by.rep.list[[i]] = curr.graph.matrix[0,]
  }
  
  for (i in 1:length(fc.heatmap.data.list)) {
    
    curr.heatmap = fc.heatmap.data.list[[i]]
    colnames(curr.heatmap) = unlist(lapply(exp.names, function(x) paste0(c("gen", x), collapse=".")))
    
    # At the end of this you'll have distributed one row from each heatmap to it's appropriate replicated
    for (j in 1:nrow(curr.heatmap)) {
      
      #curr.rep.heatmap = fc.heatmap.data.by.rep.list[[j]]
      
      #curr.rep.heatmap = rbind(curr.rep.heatmap, curr.heatmap[j,])
      
      fc.heatmap.data.by.rep.list[[j]] = rbind(fc.heatmap.data.by.rep.list[[j]], curr.heatmap[j,])
      #fc.heatmap.data.by.rep.list[[j]] = curr.rep.heatmap
      
      if (i == length(fc.heatmap.data.list)) {
        # You're on the last set so now you can name all the rows
        rownames(fc.heatmap.data.by.rep.list[[j]]) = unlist(strain.list)
        # Now you should try to plot all this data?
        
        if (j == nrow(curr.heatmap) && (j > num.rep)) {
          
          heatmap.file = paste0(c(output.file, "avg.rep.heatmap.by.fold-change"), collapse="")
          
        } else {        
          heatmap.file = paste0(c(output.file, rownames(curr.heatmap)[j], ".rep.heatmap.by.fold-change"), collapse="")
        }
        
        # heatmap.file = paste0(c(output.file, rownames(curr.heatmap)[j], "rep.heatmap.by.fold-change"), collapse=".")
        generate.strain.fold.change.heatmap(t(fc.heatmap.data.by.rep.list[[j]]), heatmap.file, "strain", "generation", "fold-change\n")
        
        curr.data = cbind(rownames(fc.heatmap.data.by.rep.list[[j]]), fc.heatmap.data.by.rep.list[[j]])
        colnames(curr.data)[1] = "strain"
        write.table(curr.data, file=paste0(c(heatmap.file, "data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
      }      
    }
  }
  
  
  
  # 6: Can we generate a fold-change comparison? What is the rate of change at each timepoint?
  # Does this normalize our data a little better?
  # Problems? If VC20019 is on a downward trajectory and the strain to compare is going up, then we'll get a large negative
  # If conversely sloping we'll get a very small negative value
  # Note: We can't do this with normalized data
  # Try to normalize data first by finding the VC20019 replicates and then applying to all the other strains?
  # fc.heatmap.data.list
  # strain.list
  
  
  # Try to get fold-change rate first
  
  # 180426: consider generating a dot-plot of this data instead? Or a line graph? Could it identify outliers?
  # The data is formatted correctly, you just need to do the stats on it at this point
  # Can you normalize the data? Inherently it should be normalized when you look at rate of change rather than total change.
  # Generate a single table of average fold-change rate with each condition as a row, each strain as a column
  
  # 180520: Overcoming the issue of when samples go only in the negative direction. Theoretically the average rate of change is -1/num.gen
  # That means if the same strain is always going downwards, the rate of change will always be the same across many samples.
  # Instead, for the average, we should stop at the first zero-valued timepoint and make that the average.
  
  
  fc.rate.data.list = list()
  
  fc.rate.avg.all = fc.heatmap.data.list[[1]][,0]
  fc.rate.avg.all = as.data.frame(cbind(rownames(fc.rate.avg.all), fc.rate.avg.all))
  colnames(fc.rate.avg.all) = c("series")

  # Generate a construct of information that makes it easier to manipulate later
  # 1) all the average data; 2) z-score; 3) z-score p-value; 4) modified z-score 5) ESD outlier analysis 
  
  fc.rate.avg.list = list()
  
  # 180606: Can we analyse based on just the final fold-change from timepoint 1 to X?
  # We can certainly do this but it may not tell us when strains drop out quickly vs slowly. We could certainly take snapshots of each 
  # expansion but I don't think this will end up being as useful.
  fc.total.list = list()
  
  for (i in 1:4) {
    fc.rate.avg.list[[i]] = fc.rate.avg.all
    fc.total.list[[i]] = fc.rate.avg.all
  }
  
  fc.rate.avg.list[[5]] = data.frame()
  fc.total.list[[5]] = data.frame()
  
  #---------- variables for Outlier analysis ----------#
  # # Minimum number of reps or samples we want to look through for outliers? I'd say 15 based on the simulation run by others
  # min.samples = 15
  # 
  # # We expect that at most, in the sample tehre will be 25% outliers
  # max.outliers = ceiling(0.25*nrow(fc.heatmap.data.list[[1]]))
  # 
  # # What should the alpha be set to for significance?
  # outlier.alpha = 0.01
  # alpha.list = c(0.1, 0.05, 0.01)
  percent.outliers = 0.25
  
  for (i in 1:length(strain.list)) {
    
    lapply(c(paste0("fold-change rate analysis: curr.strain is ", strain.list[i])), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.fc.graph = fc.heatmap.data.list[[i]]
    curr.fc.rate = curr.fc.graph[,0]
    
    fc.total.list[[1]][strain.list[i]] = as.numeric(unlist(curr.fc.graph[,ncol(curr.fc.graph)]))
    
    for (j in 2:length(exp.names)) {
      
      # 180520 we currently title our data as the fold change rate from gen.x.to.y
      # BUT we are actually dividing by the number of generations separating x and y
      # Therefore the actual calculation ends up being the average fold-change rate at each generation from x to y
      # This isn't necessarily hugely accurate when calculating our averages
      # We need to fix our definition of things when we do the average, but for now let's keep in mind that we are looking at average rate of change
      # per generation, maybe fix the x-axis on the heat map?
      
      # These numbers are already normalized to original abundance 
      
      curr.gen.diff = as.numeric(exp.names[j]) - as.numeric(exp.names[j-1])
      curr.gen.diff.title = paste0(c("gen", exp.names[j-1], "to", exp.names[j]), collapse=".")
      #curr.fc.rate[,j-1] = (curr.fc.graph[,j] - curr.fc.graph[,j-1])/curr.gen.diff
      curr.fc.rate = cbind(curr.fc.rate, ((curr.fc.graph[,j] - curr.fc.graph[,j-1])/curr.gen.diff))
      colnames(curr.fc.rate)[j-1] = curr.gen.diff.title
      
      # We need to somehow make fold-changes in the downward direction more obviously biologically relevant or equal to fold changes
      # in the positive direction. Let's alter the formula a bit to be (B-A)/A so we are always comparing to the original starting point for fold-change
      # How will this affect how the data looks? Furthermore, we will take log2 of the data so that a 10-fold vs -10-fold change look similar in -scale-
      
    }

    # 180520 replace the rowMeans calculation with a more comprehensive algorithm for defining rate of change (especially in the negative direction)
    
    # 181009: Needed to fix this average calculation as it still isn't correct in how it does it.
    # We need to send the curr.fc.graph but that's already been calculated in log2 fold-change as well so you cannot directly identify if the population has bottomed out
    # If, however, we do some mental math, the calculation of log2 fold-change is log2(curr/start) = log2(curr) - log2(start)
    # We can therefore piece together the orginal value as log2(start) + log2(curr) to see if we are below log2(min.abundance.value)
    # Fix this in get.fc.rate.avg but we need the original start population values in order to achieve this too...
    # We can't retrieve it from the strain.graph.data.list since that information is already changed to overall fold-change.
    # Solution: Pass along the original start population in log2 format. That's the only real solution to do this analysis at this point in real-time
    # Other options would be to go further back in the program and to pass a lot more information forward from the point of generating the line graphs.
    
    # PROBLEM: curr.fc.rate is already the calculated rate of change and not the log2(abundance). We need to go back further to curr.fc.rate.
    # curr.fc.graph is the answer!
    
    
    
    
    lapply(c(paste0("generating curr.fc.rate.avg")), write, output.info, append=TRUE, ncolumns=1000)
    curr.fc.rate = cbind(curr.fc.rate, get.fc.rate.avg(curr.fc.rate, curr.fc.graph, log.start.pop[i]))
    
    fc.rate.data.list[[i]] = curr.fc.rate
    
    heatmap.file = paste0(c(output.file, strain.list[i], ".rep.heatmap.by.fold-change.rate"), collapse="")
    generate.strain.fold.change.heatmap(t(fc.rate.data.list[[i]]), heatmap.file, "series", "generation difference", "rate of fold-change\n")
    
    curr.fc.rate = cbind(rownames(curr.fc.rate), curr.fc.rate)
    colnames(curr.fc.rate)[1] = "series"
    colnames(curr.fc.rate)[ncol(curr.fc.rate)] = "avg"
    write.table(curr.fc.rate, file=paste0(c(heatmap.file, "data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

    fc.rate.avg.list[[1]][strain.list[i]] = as.numeric(unlist(curr.fc.rate[,"avg"]))
    
    # generate.norm.analysis will make q-q plots of the data, zscores, modified zscore, and perhaps violin plot some of the data
    fc.rate.avg.list = generate.norm.analysis(as.numeric(curr.fc.rate[,"avg"]), "\naverage fold change", strain.list[i], fc.rate.avg.list,
                                           paste0(c(output.file, strain.list[i], ".avg.fc"), collapse=""),
                                           #alpha.list, max.outliers, min.samples, group.data)
                                           percent.outliers, group.data)
    

    # generate the norm analysis on the total/absolute fold change for strain[i]
    fc.total.list = generate.norm.analysis(as.numeric(curr.fc.graph[,ncol(curr.fc.graph)]), "\ntotal fold change", strain.list[i], fc.total.list,
                                           paste0(c(output.file, strain.list[i], ".total.fc"), collapse=""),
                                           #alpha.list, max.outliers, min.samples, group.data)
                                           percent.outliers, group.data)

  }

  # Write all of the data to a single .xlsx file
  sheet.postnames = c("data", "zscore", "zscore.pval", "zscore.mod", "outlier.info", "melt")
  # Melt that data for later generation of graphs and send it up to be stored in a single .xlsx file?
  #melt.data = list()
  #melt.data[[1]] = merge(melt(fc.rate.avg.list[[1]]), group.data, by.x="series", by.y="rep.name")
  #melt.data[[2]] = merge(melt(fc.total.list[[1]]), group.data, by.x="series", by.y="rep.name")
  fc.rate.avg.list[[6]] = merge(melt(fc.rate.avg.list[[1]]), group.data, by.x="series", by.y="rep.name")
  fc.total.list[[6]] = merge(melt(fc.total.list[[1]]), group.data, by.x="series", by.y="rep.name")
  
  
  # Now that we've generated all this aggregated normalized data analysis, let's generate some violin plots of it
  generate.norm.files(fc.rate.avg.list, output.file, "fc.rate.avg", sheet.postnames)
  generate.norm.files(fc.total.list, output.file, "fc.total", sheet.postnames)
  
  # 7: Generate heatmaps of fold-change rate data on a per-replicate basis. 
  # We need to rifle through fc.heatmap.data.list and re-assemble from X heatmaps (num strains) to Y heatmaps (num replicates)
  # One heatmap for each replicate ordered by strain (rows)
  
  fc.rate.data.by.rep.list = list()
  
  # First we have to generate place holders for each matrix in the list
  curr.graph.matrix = fc.rate.data.list[[1]]
  #colnames(curr.graph.matrix) = unlist(lapply(exp.names, function(x) paste0(c("gen", x), collapse=".")))
  #rep.names = rownames(curr.graph.matrix)
  
  
  for (i in 1:nrow(curr.graph.matrix)) {
    fc.rate.data.by.rep.list[[i]] = curr.graph.matrix[0,]
  }
  
  for (i in 1:length(fc.rate.data.list)) {
    
    curr.heatmap = fc.rate.data.list[[i]]
    #colnames(curr.heatmap) = unlist(lapply(exp.names, function(x) paste0(c("gen", x), collapse=".")))
    
    # At the end of this you'll have distributed one row from each heatmap to it's appropriate replicated
    for (j in 1:nrow(curr.heatmap)) {
      
      #curr.rep.heatmap = fc.heatmap.data.by.rep.list[[j]]
      
      #curr.rep.heatmap = rbind(curr.rep.heatmap, curr.heatmap[j,])
      
      fc.rate.data.by.rep.list[[j]] = rbind(fc.rate.data.by.rep.list[[j]], curr.heatmap[j,])
      #fc.heatmap.data.by.rep.list[[j]] = curr.rep.heatmap
      
      if (i == length(fc.rate.data.list)) {
        # You're on the last set so now you can name all the rows
        rownames(fc.rate.data.by.rep.list[[j]]) = unlist(strain.list)
        # Now you should try to plot all this data?
        
        if (j == nrow(curr.heatmap) && (j > num.rep)) {
          
          heatmap.file = paste0(c(output.file, "avg.rep.heatmap.by.fold-change.rate"), collapse="")
          
        } else {        
          heatmap.file = paste0(c(output.file, rownames(curr.heatmap)[j], ".rep.heatmap.by.fold-change.rate"), collapse="")
        }
        generate.strain.fold.change.heatmap(t(fc.rate.data.by.rep.list[[j]]), heatmap.file, "strain", "generation difference", "rate of fold-change\n")
        
        curr.data = cbind(rownames(fc.rate.data.by.rep.list[[j]]), fc.rate.data.by.rep.list[[j]])
        colnames(curr.data)[1] = "strain"
        write.table(curr.data, file=paste0(c(heatmap.file, "data.txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
        
      }      
    }
  }
  
  
}


# This function is supposed to streamline the process of choosing different metrics and analysing them
# It will generate 
# 1) q-q plots of the data
# 2) z-scores
# 3) p-values for the z-scores
# And we could Bonferroni correct it? TO do this we need to set an alpha and then divide by n (number of comparisons) to set our new p-value threshold
    # We also need to consider how to determine p-values when we have n < 30. It's suggested that we try to do a t-distribution conversion for p-values
    # 2*pt(-abs(z), df=(length(strain.list)-1))

# 4) modified z-scores that are especially useful for data that might be distrubuted with long tails
    # Publications from Iglewicz and Hoaglin suggests that we use a modified z-score instead of just a z-score

# 5) Generalized Extreme Studentized Deviate (ESD) tests for outliers
    # 180517: Run a generalized ESD test for outliers.
    # We want to do this at this point because I think we really should compare within a pool and not amongst the combined data. 
    # This is an alternative to looking at Z-scores or modified z-scores to determine outliers and it sends back a complicated object
    # Assumption for k (max outliers) would be no more than 25% of our data could be outliers, otherwise I'd really suspect the data was off
    # We can either rifle through the entire fc.rate.avg.all data frame and check each column, or we can do it as it was created. Need to decide if we can do an apply all to the function.

# 6) Violin plots of the z-score and modified z-score data that will highlight samples that are 3 (z-score) and 3.5 (modified z-score) outside of the population
generate.norm.analysis = function(dataset, data.desc, strain, data.output, output.file, percent.outliers, group.data=NULL) {
  
  # 180615 change max.outliers in function call to a percent outliers for more flexibility in using for other parts of analysis
  # Also set a default alpha.list and min.samples as these shouldn't change between different kinds of analysis
  max.outliers = ceiling(percent.outliers*length(dataset))
  alpha.list=c(0.1, 0.05, 0.01)
  min.samples = 15
  #data.output: 1) all the total fold change data; 2) z-score; 3) z-score p-value; 4) modified z-score 5) ESD outlier analysis 

  #final.output = list()
  final.output = data.output
  
  # Make a q-q plot of dataset so we can see what it looks like
  make.qq.plot(dataset, output.file, paste0(c("Normal Q-Q Plot:", strain, data.desc), collapse=" "))
  
  # Now that we've calculated the average fc.rate across the experiment we can calculate z-scores for these samples
  # Calculate the z-scores of the dataset
  final.output[[2]][strain] = scale(dataset, center=TRUE, scale=TRUE)
  
  # and the p-value of those scores
  final.output[[3]][strain] = 2*pnorm(-abs(as.numeric(unlist(final.output[[2]][strain]))))
  
  # Publications from Iglewicz and Hoaglin suggests that we use a modified z-score instead of just a z-score
  curr.mad = mad(dataset, na.rm=TRUE, high=TRUE)
  curr.median = median(dataset)
  final.output[[4]][strain] = 0.6745*(dataset - curr.median)/curr.mad
  
  # ESD test if there is data to test
  if (length(dataset) > min.samples) {
    
        for (i in 1:length(alpha.list)) {
        
        outlier.alpha = alpha.list[i]
        
        lapply(c(paste0("Attempting rosnerTest on ", strain, " with ", max.outliers, " outliers and alpha = ", outlier.alpha), collapse=""), write, output.info, append=TRUE, ncolumns=1000)
        lapply(c(paste0(dataset)), write, output.info, append=TRUE, ncolumns=1000)
        
        # ESD on the dataset provided. Try to catch any errors!
        curr.outlier.info = tryCatch(rosnerTest(dataset, k=max.outliers, alpha=outlier.alpha), error=function(e) print("Unknown failure"))

        if (curr.outlier.info == "Unknown failure") {

          lapply(c(paste0("rosnerTest failure for uknown reasons"), collapse=""), write, output.info, append=TRUE, ncolumns=1000)
          
        } else {
                  
          # Add a column to identify strain and alpha level tested
          curr.outlier.info = cbind(rep(strain, max.outliers), rep(outlier.alpha, max.outliers), curr.outlier.info$all.stats)
          
          if (!is.null(group.data)) {
            # Let's associate the group data with it as well so it's easier to identify outlier conditions
            curr.outlier.info = cbind(curr.outlier.info, group.data$rep.name[curr.outlier.info[,7]])
            # Put it into one growing data frame by concatenating the rows
    
            colnames(curr.outlier.info) = c("strain", "alpha", "rep", "rep.mean", "rep.SD", "rep.value", "obs.num", "R.i+1", "lambda.i+1", "outlier", "rep.cond")  
            
          } else {
            # There is no group.data provided and we can't do the proper assignments
            colnames(curr.outlier.info) = c("strain", "alpha", "rep", "rep.mean", "rep.SD", "rep.value", "obs.num", "R.i+1", "lambda.i+1", "outlier")  
          }
        
          final.output[[5]] = rbind(final.output[[5]], curr.outlier.info)  
        }
      }
  }

  return(final.output)
  
}

generate.norm.files =function(data.list, output.file, data.desc, sheet.info = c("data", "zscore", "zscore.pval", "zscore.mod", "outlier.info"), xlab="strain") {
  
  # Print an xlsx file of the data
  list2xlsx(data.list, paste0(c(output.file, data.desc, ".analysis.xlsx"), collapse=""), data.desc, sheet.info)

  # Print violin plots of the z-score and modified z-score data
  generate.violin(data.list[[2]], paste0(c(output.file, data.desc, ".zscore.violin.png"), collapse=""), 
                  xlab, "z-score", paste0(c("z-score violin plot of ", data.desc), collapse=""),
                  3, -3)
  generate.violin(data.list[[4]], paste0(c(output.file, data.desc, ".zscore.mod.violin.png"), collapse=""), 
                  xlab, "modified z-score", paste0(c("modified z-score violin plot of ", data.desc), collapse=""), 
                  3.5, -3.5)
  
}

# Generalized function to take in a list of input data and place it all into a .xlsx file for later analysis by the user
list2xlsx = function(input.list, output.file, data.desc, sheet.info) {

  # Write all of the data to a single .xlsx file with each list as a different sheet
  # Don't write any empty sheets
  for (i in 1:length(input.list)) {
    
    if (nrow(input.list[[i]]) > 0) {
      
      write.xlsx(input.list[[i]], file=output.file, sheetName=paste0(c(data.desc, sheet.info[i]), collapse="."), append=TRUE, row.names=FALSE)

    }
  }
}


# In order to print the violin plots, the data needs to be melted into just a couple of columns
# Right now we allow an upper and lower outlier but we should consider showing alpha level cutoffs instead?
# p=0.10 (1.65) on a z-score
# p=0.05 (1.96) on a z-score
# p=0.01 (2.58) on a z-score
generate.violin = function(data.df, plot.file, x.lab, y.lab, title, upper.outlier, lower.outlier, melted=FALSE, sigma.lines=TRUE) {
  
  if (!melted) {
    melt.data.df = melt(data.df)
    
  } else {melt.data.df = data.df}
  
  png(filename=plot.file, 
      type="cairo",
      units="in", 
      width=30, 
      height=10, 
      pointsize=font.size, 
      res=300   )

  # vplot = ggplot(melt.data.df[melt.data.df$value > lower.outlier & melt.data.df$value < upper.outlier,], aes(x=variable, y = value)) + 
  #   xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
  #   geom_violin() + theme(text = element_text(size=12), axis.text.x = element_text(angle = 90, hjust=1)) + 
  #   (geom_jitter(data=melt.data.df[melt.data.df$value <=lower.outlier | melt.data.df$value >= upper.outlier,],shape=16, position = position_jitter(0.5), aes(colour="red"))) + 
  #   geom_jitter(data=melt.data.df[melt.data.df$value >lower.outlier & melt.data.df$value < upper.outlier,],shape=16, position = position_jitter(0.5)) +
  #   theme(legend.position="bottom") + scale_colour_discrete(name="", labels="Candidate\nOutliers")

  if (nrow(melt.data.df[melt.data.df$value <=lower.outlier | melt.data.df$value >= upper.outlier,]) > 0) {

    vplot = ggplot(melt.data.df[melt.data.df$value > lower.outlier & melt.data.df$value < upper.outlier,], aes(x=variable, y = value)) + 
      xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
      geom_violin() + theme(text = element_text(size=12), axis.text.x = element_text(angle = 90, hjust=1)) + 
      (geom_jitter(data=melt.data.df[melt.data.df$value <=lower.outlier | melt.data.df$value >= upper.outlier,],shape=16, position = position_jitter(0.5), aes(colour="red"))) + 
      geom_jitter(data=melt.data.df[melt.data.df$value >lower.outlier & melt.data.df$value < upper.outlier,],shape=16, position = position_jitter(0.5)) +
      theme(legend.position="bottom") + scale_colour_discrete(name="", labels="Candidate\nOutliers")
  } else {
    # There are no outliers and adding them with different colours will break ggplot
    vplot = ggplot(melt.data.df[melt.data.df$value > lower.outlier & melt.data.df$value < upper.outlier,], aes(x=variable, y = value)) + 
      xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
      geom_violin() + theme(text = element_text(size=12), axis.text.x = element_text(angle = 90, hjust=1)) + 
      geom_jitter(data=melt.data.df[melt.data.df$value >lower.outlier & melt.data.df$value < upper.outlier,],shape=16, position = position_jitter(0.5))

  }

  if (sigma.lines) {
    vplot = vplot + 
      geom_hline(yintercept=1.65, linetype="dashed", color="yellow") + geom_hline(yintercept=-1.65, linetype="dashed", color="yellow") +
      geom_hline(yintercept=1.96, linetype="dashed", color="orange") + geom_hline(yintercept=-1.96, linetype="dashed", color="orange") +
      geom_hline(yintercept=2.58, linetype="dashed", color="red") + geom_hline(yintercept=-2.58, linetype="dashed", color="red") 
  }
  print(vplot)
  dev.off()

  
}

### ----------------------- Correlation analysis of data -----------------------###

# Generates a correlation of the desired type (pearson or spearman)
# Right now it will default to pearson as this is the most appropriate looking
correlate.data = function(data, cor.type) {
  
  cor.data = cor(data, method=cor.type)
  #cor.data = cor(data, method="pearson")
  
  return (cor.data)
}

# Given a data matrix of strains (rows) and replicates (cols)
# Note that data.matrix must be in the matrix format.
generate.norm.heatmap = function(data, norm.strain, output.file) {

  graph.title = tail(unlist(strsplit(output.file, "\\/")), 1)
  
  # normalize the data and it will return a matrix
  norm.matrix = normalize.strains(data, norm.strain)
  
  ramp = colorRamp(c("white", "red","yellow","orange","green", "blue"))
  
  # values/settings
  # at=c(a,b,c) allows us to set a range or series of ranges.
  # at=c(seq(0,10,0.01),seq(10.01,30,19.99)) lets us set a range from 0-10, and 10.01-30. 
  # The last range section ends up being coloured similarly allowing more colour gradient in the 0-10 range
  # We should scale the output for larger files if we have more replicates to show. (used to be 16.5w x 15h) 
  # We determine width by col and height by rows since we'll be using the transpose of the matrix to plot.
  
  #img.width = max(ncol(norm.matrix)*font.size*3/300, 15)
  #img.height = max(nrow(norm.matrix)*font.size*3/300, 15)
  
  img.width = max((ncol(norm.matrix)*50 + 300)/300, 4)
  img.height = max((nrow(norm.matrix)*50 + 300)/300, 4)
  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.height, 
      pointsize=font.size, 
      res=300   )
  
  data.plot = levelplot(t(norm.matrix), main=graph.title, xlab="series", ylab="strain", 
                        col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                        cuts=10000, at=c(seq(0,10,0.01),seq(10.01,30,19.99)), scales=list(x=list(rot=90)))
  print(data.plot)
  dev.off()
  
  return (list(data.plot, norm.matrix))
  
}

# Given a data matrix of strains (rows) and replicates (cols)
# Note that data.matrix must be in the matrix format.
generate.strain.heatmap = function(data, norm.strain=NULL, output.file) {
  
  output.file.norm = paste(c(output.file, ".norm.png"), collapse="")
  graph.title.norm = tail(unlist(strsplit(output.file.norm, "\\/")), 1)
  
  output.file = paste(c(output.file, ".png"), collapse="")  
  graph.title = tail(unlist(strsplit(output.file, "\\/")), 1)
  
  data = as.matrix(data)
  
  # values/settings
  # at=c(a,b,c) allows us to set a range or series of ranges.
  # at=c(seq(0,10,0.01),seq(10.01,30,19.99)) lets us set a range from 0-10, and 10.01-30. 
  # The last range section ends up being coloured similarly allowing more colour gradient in the 0-10 range
  # We determine width by col and height by rows since we'll be using the transpose of the matrix to plot.
  
  ramp = colorRamp(c("white", "red", "orange","yellow", "green", "blue"))
  #ramp = colorRamp(c("black", "white", "red","yellow","orange","green", "blue"))
  #img.width = max(ncol(data)*font.size*3/300, 15)
  #img.height = max(nrow(data)*font.size*3/300, 15)
  
  img.width = max((ncol(data)*50 + 300)/300, 4)
  img.height = max((nrow(data)*50 + 300)/300, 4)
  
  if (class(norm.strain) != "NULL") {
    # normalize the data and it will return a matrix
    norm.matrix = normalize.strains(data, norm.strain)
     
    png(filename=output.file.norm, 
        type="cairo",
        units="in", 
        width=img.width, 
        height=img.height, 
        pointsize=font.size, 
        res=300   )
    
    # Round to the nearest tens digit and then go upwards if necessary
    # We should also check the case where our VC20019 has gone to 0 for whatever reason. This will generate an "Inf" value we'd like to get rid of before
    max.range = round(max(norm.matrix[(norm.matrix < Inf)], na.rm = TRUE), digits = -1)
    
    
    # Special case on how to make this if the max.range is less than 10
    if (max.range < 10) {
      
      max.range = round(max(norm.matrix[(norm.matrix < Inf)], na.rm = TRUE))
      
      if (max.range < max(norm.matrix[(norm.matrix < Inf)], na.rm = TRUE)) {
        max.range = max.range + 1
      }
      
    } else if (max.range < max(norm.matrix[(norm.matrix < Inf)], na.rm = TRUE)) {
      
      max.range = max.range + 10
    }
    
    if (max.range <= 10) {
      data.plot.norm = levelplot(t(norm.matrix), main=graph.title.norm, xlab="series", ylab="strain", ylab.right=paste0(c("Normalized abundance vs. Control (", norm.strain,")\n"), collapse=""), 
                                 col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                                 cuts=10000, at=c(seq(0,max.range,0.01)), scales=list(x=list(rot=90)),
                                 colorkey = list(at=c(seq(0,max.range,0.01)), labels=list(at=c(seq(0,max.range,1)), labels=c(seq(0,max.range,1))))
      )
      print(data.plot.norm)
    } else {
      
      # I've played around with this command call in order to fix the labels and coloring
      # the colorkey command mirrors the "at" commend above it but allows me to label it so that above 10, we don't have to show how high it stretches
      data.plot.norm = levelplot(t(norm.matrix), main=graph.title.norm, xlab="series", ylab="strain", ylab.right=paste0(c("Normalized abundance vs. Control (", norm.strain,")\n"), collapse=""), 
                                 col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                                 cuts=10000, at=c(seq(0,10,0.01),seq(10.01,max.range,max.range-10.01)), scales=list(x=list(rot=90)),
                                 colorkey = list(at=c(seq(0,10,0.01),seq(10.01,12,1.99)), labels=list(at=c(seq(0,8,2), 10, 12), labels=c(seq(0,8,2), ">10", max.range)))
      )
      print(data.plot.norm)
    }
    
    dev.off()

  } else {
    data.plot.norm=NULL
    norm.matrix=NULL
  }
  # Generate the raw abundance version. This should have no NA values in it in theory?
  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.height, 
      pointsize=font.size, 
      res=300   )
  
  #max.range = 1
  max.range = round(max(data, na.rm = TRUE), digits = 1)
  if (max.range < max(data, na.rm = TRUE)) {
    
    max.range = max.range + 0.1
  }
  
  data.plot = levelplot(t(data), main=graph.title, xlab="series", ylab="strain", ylab.right="abundance\n",
                        col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                        cuts=10000, at=c(seq(0,max.range,0.001)), scales=list(x=list(rot=90)))
                        #cuts=10000, at=c(seq(0,max.range,0.01)), scales=list(x=list(rot=90)))
                        #cuts=10000, at=c(seq(0,0.01,0.01), seq(0.011,max.range,0.01)), scales=list(x=list(rot=90)))
  print(data.plot)
  dev.off()
  
  
  return (list(data.plot.norm, data.plot, norm.matrix))
  
}

generate.strain.fold.change.heatmap = function(data, output.file, x.title, y.title, y.legend) {
  
  output.file = paste(c(output.file, ".png"), collapse="")  
  graph.title = tail(unlist(strsplit(output.file, "\\/")), 1)
  
  data = data.matrix(data)
  
  # values/settings
  # at=c(a,b,c) allows us to set a range or series of ranges.
  # at=c(seq(0,10,0.01),seq(10.01,30,19.99)) lets us set a range from 0-10, and 10.01-30. 
  # The last range section ends up being coloured similarly allowing more colour gradient in the 0-10 range
  # We determine width by col and height by rows since we'll be using the transpose of the matrix to plot.
  
  ramp = colorRamp(c("white", "red", "orange","yellow", "green", "blue"))
  #ramp = colorRamp(c("black", "white", "red","yellow","orange","green", "blue"))
  #img.width = max(ncol(data)*font.size*3/300, 15)
  #img.height = max(nrow(data)*font.size*3/300, 15)
  
  img.width = max((ncol(data)*50 + 300)/300, 4)
  img.height = max((nrow(data)*50 + 300)/300, 4)
  
  # Generate the raw abundance version. This should have no NA values in it in theory?
  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.height, 
      pointsize=font.size, 
      res=300   )
  
  max.range = round(max(data, na.rm = TRUE), digits = 1)
  if (max.range < max(data, na.rm = TRUE)) {
    
    # This could prematurely set the range to 1. When this is called, what are the number ranges?
    # Change how we'll manipulate it depending on if the scale is 0-1 or 1-X
    if (max.range < 1) {
      max.range=max.range + 0.1
    } else {
      max.range = max.range + 1  
    }
    
  }
  
  min.range = round(min(data, na.rm = TRUE), digits = 1)
  if (min.range >= 0) {min.range = 0} 
  else if (min.range > min(data, na.rm = TRUE)) {
    
    # Does this ever happen? Not using current calculations for the heatmap. Maybe if fold-change is calculated differently...
    min.range = min.range - 1
  }
  
  data.plot = levelplot(t(data), main=graph.title, xlab=x.title, ylab=y.title, ylab.right=y.legend,
                        col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                        cuts=10000, at=c(seq(min.range,max.range,((max.range-min.range)/10000))), scales=list(x=list(rot=90)))
  #cuts=10000, at=c(seq(0,max.range,0.01)), scales=list(x=list(rot=90)))
  #cuts=10000, at=c(seq(0,0.01,0.01), seq(0.011,max.range,0.01)), scales=list(x=list(rot=90)))
  print(data.plot)
  dev.off()
  
  
  return (list(data.plot))
  
}

# Given a data set of strains (row) x replicates (col), generate a correlation of the pools
# Note that data is expected to be in a data frame
# If the data frame has NA values or INF values how are they plotted?
generate.cor.heatmap = function(pool.data.df, cor.method, output.file) {
  
  output.file = paste0(c(output.file, ".png"), collapse="")
  graph.title = tail(unlist(strsplit(output.file, "\\/")), 1)
  
  # normalize the data and it will return a matrix
  cor.matrix = correlate.data(pool.data.df, cor.method)
  
  #ramp = colorRamp(c("white", "red","yellow","orange","green", "blue"))
  #ramp = colorRamp(c("blue", "green","orange","yellow", "red"))
  ramp = colorRamp(c("blue", "green", "red"))
  # values/settings
  # at=c(a,b,c) allows us to set a range or series of ranges.
  # at=c(seq(0,10,0.01),seq(10.01,30,19.99)) lets us set a range from 0-10, and 10.01-30. 
  # The last range section ends up being coloured similarly allowing more colour gradient in the 0-10 range
  
  # Round down to the lowest tenth
  min.range = round(min(cor.matrix, na.rm = TRUE), digits = 1)
  #lapply(c("Generating cor heatmap with min(cor.matrix): ", min(cor.matrix)), write, output.info, append=TRUE, ncolumns=1000)
  #lapply(c("Generating cor heatmap with min.range: ", min.range), write, output.info, append=TRUE, ncolumns=1000)
  if (min.range > min(cor.matrix, na.rm = TRUE)) {
    
    min.range = min.range - 0.1
  }
  
  # Replace the values of 22x22 with adjusted values
  # img.width = max(ncol(cor.matrix)*font.size*3/300, 22)
  
  img.width = max((ncol(cor.matrix)*50 + 300)/300, 15)

  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.width, 
      pointsize=font.size, 
      res=300   )
  
  data.plot = levelplot(t(cor.matrix), main=graph.title, xlab="series", ylab="series", ylab.right="Correlation Coefficient\n",
                        col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                        cuts=10000, at=c(seq(min.range,1,0.01), 1.01, Inf), scales=list(x=list(rot=90)))
  print(data.plot)
  dev.off()
  
  # 181126: Can we make a dendrogram of this data for analysis?
  # We can't just do a dendrogram on the spearman data but should rather do it based on all of the abundance data within a sample? 
  # Looking at M10 and M11, the spearman data seems to line up better as groups than just giving it the raw data but is that correct?
  # 181127: Spoke with Lina and she suggested it was OK to do the Spearman before feeding to heatmap.2
  # Can play with the dendrogram and give it number of groups too but altering hclust function
  
  output.file = paste0(c(output.file, "dendrogram.png"), collapse="")
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.width, 
      pointsize=font.size, 
      res=300   )
  
  #heatmap.2(cor.matrix, trace="none", dendrogram="column")
  plot(hclust(as.dist(1-t(cor.matrix)), method="complete"))
  dev.off()
  
  
  
  return (list(data.plot, cor.matrix))
  
}

# Generic heatmap function
generate.generic.heatmap = function(data.df, output.file) {
  
  output.file = paste0(c(output.file, ".png"), collapse="")
  graph.title = tail(unlist(strsplit(output.file, "\\/")), 1)
  
  data.matrix = as.matrix(data.df)
  
  #ramp = colorRamp(c("white", "red","yellow","orange","green", "blue"))
  #ramp = colorRamp(c("blue", "green","orange","yellow", "red"))
  ramp = colorRamp(c("blue", "green", "red"))
  # values/settings
  # at=c(a,b,c) allows us to set a range or series of ranges.
  # at=c(seq(0,10,0.01),seq(10.01,30,19.99)) lets us set a range from 0-10, and 10.01-30. 
  # The last range section ends up being coloured similarly allowing more colour gradient in the 0-10 range
  
  # Round down to the lowest tenth
  min.range = round(min(data.matrix, na.rm = TRUE), digits = 1)
  max.range = round(max(data.matrix, na.rm = TRUE), digits = 1)
  #lapply(c("Generating cor heatmap with min(cor.matrix): ", min(cor.matrix)), write, output.info, append=TRUE, ncolumns=1000)
  #lapply(c("Generating cor heatmap with min.range: ", min.range), write, output.info, append=TRUE, ncolumns=1000)
  if (min.range > min(data.matrix, na.rm = TRUE)) {
    
    min.range = min.range - 0.1
  }
  if (max.range < max(data.matrix, na.rm = TRUE)) {
    max.range = max.range + 0.1
  } 
  
  # Replace the values of 22x22 with adjusted values
  # img.width = max(ncol(cor.matrix)*font.size*3/300, 22)
  
  img.width = max((ncol(data.matrix)*50 + 300)/300, 15)
  
  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.width, 
      pointsize=font.size, 
      res=300   )
  
  data.plot = levelplot(t(data.matrix), main=graph.title, xlab="series", ylab="series", ylab.right="Correlation Coefficient\n",
                        col.regions=rgb(ramp(seq(0,1,length=10000)), max=255), 
                        cuts=10000, at=c(seq(min.range,max.range,0.01), (max.range + 0.01), Inf), scales=list(x=list(rot=90)))
  print(data.plot)
  dev.off()
  
  return (list(data.plot, data.matrix))
  
}


# 180520: This is a more comprehensive way to calculate average rate of change from a table of rate of change
# input: take in a table of rate of change, where each row is an experiment and each column is the rate of change (calculated per generation) between adjacent rows
# curr.fc.rate is a table which also has a first column of information
# Note that at current the rate of change is determined per generation and not per time point. If that is the case then
# we can divide by the number of rates of change calculated to get our average (Rather then dividing by total generations in the time course) 
# return: a vector of averages where rate of change is calculated from only the first zero entry and backwards

get.fc.rate.avg = function(curr.fc.rate, curr.fc.graph, log.start.pop) {
  
  curr.fc.rate.avg = list()
  curr.fc.ab = curr.fc.graph + log.start.pop
  
  for (i in 1:nrow(curr.fc.rate)) {
    
    # grab the current row
    curr.row = curr.fc.rate[i,]
    # regenerate the original abundance values
    curr.row.log.abundance = curr.fc.ab[i,]
    
    # Set the default divisor for the calculation
    num.fc.rates = length(curr.row)
    
    # Check ifyou need to change the divisor by finding the first zero
    # This is actually a weird point to check since it could be a negative number if we are doing log transformations.
    # We should set a minimum value instead
    # 181009 - this is not the correct place to make this calculation.
    # The values are already converted to a log2 value that's the (B-A)/num.gens
    # What we must accomplish here is looking for 0 to no change
    # recalculate the original value! log2(start) + log2(current)

    lapply(c(paste0("curr.row: ")), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0(curr.row)), write, output.info, append=TRUE, ncolumns=1000)

    
    lapply(c(paste0("curr.row.log.abundance: ")), write, output.info, append=TRUE, ncolumns=1000)
    lapply(c(paste0(curr.row.log.abundance)), write, output.info, append=TRUE, ncolumns=1000)

    
    if (any(curr.row.log.abundance <= min.avg.rate)) {
     
      lapply(c(paste0("Low log abundance found!")), write, output.info, append=TRUE, ncolumns=1000)
       
      # Note the first occurrence of 0. We can use this column number as the divisor
      # num.fc.rates = which(curr.row.log.abundance <= min.avg.rate)[1]
      
      # subtract 1 from the column because of how abundance is represented (a,b,c) vs how fc.rate is represented (a-b, b-c) 
      num.fc.rates = which(curr.row.log.abundance <= min.avg.rate)[1] - 1
      lapply(c(paste0("num.fc.rates is: ", num.fc.rates)), write, output.info, append=TRUE, ncolumns=1000)

    }
    
    curr.fc.rate.avg = c(curr.fc.rate.avg, sum(as.numeric(curr.row[1:num.fc.rates]))/num.fc.rates)
    
  }
  
  return(as.numeric(curr.fc.rate.avg))
  
}

make.qq.plot = function(data, output.file, title) {
  
  output.file = paste0(c(output.file, "q-q.plot.png"), collapse=".")
  img.width = 8
  
  png(filename=output.file, 
      type="cairo",
      units="in", 
      width=img.width, 
      height=img.width, 
      pointsize=font.size, 
      res=300   )
      
  data.plot = qqnorm(data, main=title)
  
  #print(data.plot)
  dev.off()

}

make.multiple.boxplot = function(data.df, file.info) {
  
  strains = unlist(unique(data.df$variable))
  
  for (i in 1:length(strains)) {
    
    file.name = paste0(c(file.info, strains[i], "png"), collapse=".")
    title = paste0(c(strains[i], "mean fold-change rate by temperature"), collapse=" ")
    
    curr.data = data.df[data.df$variable == strains[i],]
    
    make.boxplot.by.group(curr.data, file.name, title)
    
  }
  
  
}

make.boxplot.by.group = function(data.df, file.name, title) {
# 
#   vplot = ggplot(data.df[data.df$group %in% c("HT115", "OP50", "Na22"),], aes(x=group, y = value)) + 
#     xlab("Food Source") + ylab(expression('Mean fold-change rate (log'[2]*')')) + ggtitle(title) +
#     geom_boxplot(outlier.shape = NA)+theme(text = element_text(size=18), axis.text.x = element_text(angle = 0)) + 
#     scale_x_discrete(limits=c("HT115", "Na22", "OP50"), labels=c("HT115", "Na22", "OP50")) +
#     geom_jitter(data=data.df[(data.df$group == "HT115"),],shape=16, size=3, position = position_jitter(0.5)) +
#     geom_jitter(data=data.df[(data.df$group == "Na22"),],shape=16, size=3, position = position_jitter(0.5)) + 
#     geom_jitter(data=data.df[(data.df$group == "OP50"),],shape=16, size=3, position = position_jitter(0.5)) +
#     scale_color_manual(name="Food Source", values=jitter.palette, labels=c("HT115", "Na22", "OP50")) + theme(legend.position="bottom")

  vplot = ggplot(data.df[data.df$group %in% c("HT115", "Low-Temp", "High-Temp"),], aes(x=group, y = value)) + 
    xlab("Temperature") + ylab(expression('Mean fold-change rate (log'[2]*')')) + ggtitle(title) +
    geom_boxplot(outlier.shape = NA)+theme(text = element_text(size=18), axis.text.x = element_text(angle = 0)) + 
    scale_x_discrete(limits=c("Low-Temp", "HT115", "High-Temp"), labels=c("16C", "20C", "25C")) +
    geom_jitter(data=data.df[(data.df$group == "HT115"),],shape=16, size=3, position = position_jitter(0.5)) +
    geom_jitter(data=data.df[(data.df$group == "Low-Temp"),],shape=16, size=3, position = position_jitter(0.5)) + 
    geom_jitter(data=data.df[(data.df$group == "High-Temp"),],shape=16, size=3, position = position_jitter(0.5)) +
    scale_color_manual(name="Temperature", values=jitter.palette, labels=c("HT115", "Na22", "OP50")) + theme(legend.position="bottom")

  
  
  vplot2 = vplot +
      geom_hline(yintercept=-0.62, linetype="dashed", color="red", size=1) +
      geom_hline(yintercept=-0.152, linetype="dashed", color="orange", size=1) +
      geom_hline(yintercept=0.2, linetype="dashed", color="green", size=1)

  png(filename=file.name,
      type="cairo",
      units="in",
      width=8,
      height=8.5,
      pointsize=font.size,
      res=300   )

  print(vplot2)
  dev.off()
  
}


# make.boxplot.by.group = function(data.df, group.list, x.lab, y.lab, title, x.data, y.data, legend.title, legend.labels, output.file, var.list) {
#   
#   jitter.palette = c("#000000", "#FF0033", "#3399FF", "#990099", "#33FF00", "#CC9900", "#FF9966")
#   
#   for (i in 1:length(var.list)) {
#     
#     file.name = paste0(c(output.file, var.list[i], "png"), collapse=".")
#     curr.title = paste0(c(var.list[i], title), collapse=" ")
#     
#     curr.data = data.df[data.df$variable == var.list[i],]
#     
#     vplot = ggplot(curr.data[curr.data$group %in% group.list,], aes(x=group, y = variable)) + 
#       xlab(x.lab) + ylab(y.lab) + ggtitle(curr.title) +
#       geom_boxplot(outlier.shape = NA)+theme(text = element_text(size=18), axis.text.x = element_text(angle = 0, hjust=1)) +
#       scale_x_discrete(limits=group.list, labels=legend.labels)+
#       geom_jitter(data=curr.data[(curr.data$group == group.list[1]),],shape=16, size=3, position = position_jitter(0.5), aes(color=group)) +
#       geom_jitter(data=curr.data[(curr.data$group == group.list[2]),],shape=16, size=3, position = position_jitter(0.5), aes(color=group)) +
#       geom_jitter(data=curr.data[(curr.data$group == group.list[3]),],shape=16, size=3, position = position_jitter(0.5), aes(color=group))
# 
#     # 
#     # for (j in 1:length(group.list)) {
#     # 
#     #   vplot = vplot + geom_jitter(data=curr.data[(curr.data$group == group.list[j]),], shape=16+j, size=3, position = position_jitter(0.5), aes(color=group))
#     # 
#     # }
# 
#     vplot = vplot + scale_color_manual(name=legend.title, values=jitter.palette, labels=legend.labels) + theme(legend.position="bottom")
# 
#     vplot2 = vplot +
#       geom_hline(yintercept=-0.62, linetype="dashed", color="red", size=1) +
#       geom_hline(yintercept=-0.152, linetype="dashed", color="orange", size=1) +
#       geom_hline(yintercept=0.2, linetype="dashed", color="green", size=1)
# 
#     png(filename=file.name, 
#         type="cairo",
#         units="in", 
#         width=max(length(group.list)/2, 8), 
#         height=8.5, 
#         pointsize=font.size, 
#         res=300   )
#     
#     print(vplot2)
#     dev.off()
#     
#   }
#   
# }

### ----------------------- Derelict functions that have been replaced with more appropriate versions -----------------------###

# The purpose of this function is to generate "growth" curves given a series of boxplot data
# This essentially helps us to determine of a strain is growing or not
# This generates a line graph similar to the make.boxplot.growth.curve function except it really doesn't look as good. 
make.growth.curve = function(boxplot.list, strain.list, group.list, output.directory) {
  
  
  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  # Create the output directory
  dir.create(output.directory)
  
  
  quartiles.list = list()
  
  # First format the boxplot data into what you want
  # Keep them in a chained list
  for (i in 1:(length(boxplot.list))) {
    
    curr.boxplot = boxplot.list[[i]]
    
    # What if we're looking at a single expansion, no replicates? There'll be no boxplot information
    # The value is just the average obtained from a single set of MIPs????
    # If there is boxplot data, it has more than 1 row. We'll rename the columns on the data at this point
    # otherwise it should already be named (and likely is a data frame)
    if (!is.data.frame(curr.boxplot)) {
      
      curr.boxplot = curr.boxplot$stats
      colnames(curr.boxplot) = boxplot.list[[i]]$names        
    }
    
    quartiles.list[[i]] = curr.boxplot
    
  }
  
  growth.curve.list = vector("list", length(strain.list))
  
  for (i in 1:length(quartiles.list)) {
    
    curr.quartile = as.data.frame(quartiles.list[[i]])
    
    for (j in 1:length(strain.list)) {
      
      curr.strain = strain.list[j]
      growth.curve.list[[j]] = cbind(growth.curve.list[[j]], curr.quartile[,curr.strain])
      
    }
    
  }
  
  # Need to iterate through each and make it a data frame
  # Then label the columns
  # Then add a new column identifying each row as a specific expansion
  # melt the data together based on the expansion variable
  # ggplot(VC20152.melt, aes(x=expansion, y=value, group=variable, colour=variable)) + geom_line()
  
  for (i in 1:length(growth.curve.list)) {
    
    curr.file.name = paste(c(output.directory, strain.list[i], ".growth.curve.png"), collapse = "")
    #curr.file.name = paste(c(strain.list[i], "growth.curve.png"), collapse = "")
    
    curr.quart.df = as.data.frame(t(as.matrix(growth.curve.list[[i]])))
    curr.quart.df = cbind(group.list, curr.quart.df)
    colnames(curr.quart.df) = c("expansion", "lower.whisker", "lower.quart", "mean", "upper.quart", "upper.whisker")
    curr.quart.melt = melt(curr.quart.df, id.vars="expansion")
    colnames(curr.quart.melt) = c("expansion", "variable", "value")
    
    png(filename=curr.file.name, 
        type="cairo",
        units="in", 
        width=11, 
        height=8.5, 
        pointsize=14, 
        res=300   )
    
    curr.plot = ggplot(curr.quart.melt, (aes(x=expansion, y=value, group=variable, colour=variable))) + geom_line() + ggtitle(strain.list[i])
    ggsave(curr.plot, filename=curr.file.name)
    
    dev.off()    
  }
  
  
  return (growth.curve.list)
  
}

