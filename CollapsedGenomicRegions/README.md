### A workflow and Python script for finding likely collapsed regions in a genome assembly

###### Preparation
What is needed to run this workflow:
* A genome assembly FASTA file
* The long reads that were used to produce the assembly (PacBio/Nanopore)

###### Map long reads back to the assembly 
For PacBio reads use:
`minimap2 -t 4 -ax map-pb ${genome} ${longreads} --secondary=no -o mapping.sam`
For Nanopore reads use:
`minimap2 -t 4 -ax map-ont ${genome} ${longreads} --secondary=no -o mapping.sam`

For further information about minimap2 see here: https://github.com/lh3/minimap2
