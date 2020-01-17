PhenoMIPs analysis scripts and tools

Order of operations

1) Initial sequencing data should be run on PhenoMIP.analysis.V1.R. Provided are batch.run.example.txt, example.MIP.info, example.fastq.info.txt. These will be needed to generate the initial abundance calculations for all MIPs and strains.

2) After the abundance data and files are generated, they must be reformated to work with the secondary script "PhenoMIP.graph.analysis.v1.0.R". This secondary script will take in the "strain.average.txt" file as the data.file for the generate.plots() function. This function can be batch called when working with multiple data sets or different pools of strains. This step relies on a specific data format so please look at the supplied example "M1.140508.organized.renamed.txt" file for more information. The batch function example is named "example.dataset.info.txt" and each data set will also need to reference a "group" file for plotting purposes (see "M1.group.info.txt" as an example).

