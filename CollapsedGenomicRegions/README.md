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
