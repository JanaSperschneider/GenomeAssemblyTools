### A workflow and Python script for finding likely collapsed regions in a genome assembly

##### Preparation
What is needed to run this workflow:
* A genome assembly FASTA file
* The long reads that were used to produce the assembly (PacBio/Nanopore)

##### 1) Map long reads back to the assembly 
For PacBio reads use:
`minimap2 -t 4 -ax map-pb ${genome} ${longreads} --secondary=no -o mapping.sam`

For Nanopore reads use:
`minimap2 -t 4 -ax map-ont ${genome} ${longreads} --secondary=no -o mapping.sam`

For further information about minimap2 see here: https://github.com/lh3/minimap2

##### 2) Map long reads back to the assembly 
Use BBMap's pileup.sh tool to calculate binned read coverage per location, binsize=1000 should work well for most genomes.

`pileup.sh in=mapping.sam out=contig_coverage.txt basecov=basecoverage.txt binsize=1000 bincov=bincoverage.txt`

For further information about bbmap see here: https://sourceforge.net/projects/bbmap/

##### 3) Plot the coverage and decide on the collapsed region threshold
Use the supplied R script (Coverage.R) to plot the file bincoverage.txt and inspect the coverage plot.

For example, see this plot here:

<a href="url"><img src="https://github.com/JanaSperschneider/GenomeAssemblyTools/blob/master/CollapsedGenomicRegions/Coverage_Example.png" align="left" height="48" width="48" ></a>

##### 4) Run the Python script to get likely collapsed regions in your genome
