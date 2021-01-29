### A workflow for finding contaminants in a fungal genome assembly

##### Preparation
What is needed to run this workflow:
* A genome assembly FASTA file
* The long reads that were used to produce the assembly (PacBio/Nanopore) for coverage analysis
* BLAST databases

##### 1) Map long reads back to the assembly to get coverage for each contig
For PacBio (non-Hifi) reads use:
`minimap2 -ax map-pb ${genome} ${longreads} --secondary=no -o mapping.sam`

For PacBio (Hifi) reads use:
`minimap2 -ax asm20 ${genome} ${longreads} --secondary=no -o mapping.sam`

For Nanopore reads use:
`minimap2 -ax map-ont ${genome} ${longreads} --secondary=no -o mapping.sam`

For further information about minimap2 see here: https://github.com/lh3/minimap2

Use BBMap's pileup.sh tool to calculate read coverage and GC content per contig.

`pileup.sh in=mapping.sam out=contig_coverage.txt`
For further information about bbmap see here: https://sourceforge.net/projects/bbmap/

###### From the contig coverage table, one can now identify very low-coverage contigs (e.g. < 2x coverage) and remove those from the assembly.

##### 2) Identify mitochondrial contig

From the contig coverage table, the mitochondrial contig should be high coverage and also low in GC content.
Confirm this with a BLAST search against a mitochondrial database from NCBI.

`makeblastdb -in mitochondrion.1.1.genomic.fna -dbtype nucl -out MITO
blastn -query $genome -db MITO -dust yes -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | gawk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > mito.screen.txt
cat mito.screen.txt | grep -v "^#" >> MITO_CONTIGS.txt`
