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

**From the contig coverage table, one can now identify very low-coverage contigs (e.g. < 2x coverage) and remove those from the assembly.**

##### 2) Identify mitochondrial contig

From the contig coverage table, the mitochondrial contig should be high coverage and also low in GC content.
Confirm this suspect contig as mitochondrial with a BLAST search against a mitochondrial database from NCBI (https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/).

`makeblastdb -in mitochondrion.1.1.genomic.fna -dbtype nucl -out MITO`

`blastn -query $genome -db MITO -dust yes -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | gawk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' > mito.screen.txt`

`cat mito.screen.txt | grep -v "^#" >> MITO_CONTIGS.txt`

`cut -f1 MITO_CONTIGS.txt | sort | uniq -c | sort -nr | head`

**Remove the mitochondrial contig from the assembly and put into separate FASTA file.**

##### 3) Identify contaminant contigs

Run a BLAST search for all your contigs and record the hits. 

```o="_contig.screen.txt"
for f in *.fasta
	do
		n=$(basename $f)
		echo $n
		FILE=$n$o
		if [ -f "$FILE" ]; then
		    echo "$FILE exists."
		else 
		    echo "$FILE does not exist."
		    blastn -query $f -db nt -evalue 1e-5 -perc_identity 75 -outfmt "6 -qseqid sseqid pident length evalue bitscore sgi sacc staxids sscinames scomnames stitle" > $FILE
		fi	
	done
```

A quick scan can identify if the top hits map to bacteria. For example, this one-liner would output contigs in a rust fungal assembly that are likely to be contaminants:

```grep -E -m1 -v 'Puccinia|Phakopsora|psidii|Melamps|Uromyces|ribosomal' *_contig.screen.txt```

**Remove the contaminant contigs from the assembly.**




