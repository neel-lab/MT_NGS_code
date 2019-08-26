indel_parse_whole.awk: Indel Parsing AWK Script

Script Explanation:

Using the CIGAR strings from BAM files, and the chromosomal loci of the expected CRISPR on-target mutations, we deduced whether or not a read harbors a legitimate CRISPR mutation.

This code processess bed files that contain intersection information between the reads from our experiment and the CRISPR target sites.  In our manuscript, we first found all reads that intersected with CRISPR on-target sites using bedtools.  We used the bedtools "intersect" command as follows to obtain these bed files:

bedtools intersect -a bed_formatted_BAM_file.bed -b On_target_CRISPR_Sites.bed -wo > On_target_reads.bed

The "On_target_CRISPR_Sites.bed" file has the following fields:

chr start_target end_target start_sgRNA end_sgRNA strand GeneName

start_target and end_target refer to the bases 15-17 bases downstream from the start of the protospacer sequence.  This is where the indel events are expected to take place.  start_sgRNA and end_sgRNA refer to the start and end of the entire 20nt protospacer sequence.  Make sure that "On_target_CRISPR_Sites.bed" has chromosome loci referring to the last 15 to 17 bases of the protospacer sequence.  

The fields created from this bedtools intersection were as follows:

chr_read	start_read	end_read	CIGAR gene_name	read_seq	ref_seq	chr_targetSite	start_targetSite	end_targetSite	read_type	

Field Descriptions:
1. chr_read: chromosome of the read
2. start_read: read's start location on the chromosome
3. end_read: read's end location of the chromosome
4. CIGAR: alignment information of the read to the genome
5. gene_name: gene symbolic name obtained from intersecting our raw reads with a bed file of the human genome GRCh37
6. read_seq: sequence directly from the sequencer
7. ref_seq: sequence parsed from the genome
8. chr_targetSite: chromosome of the CRISPR 3bp target site (same as chr_read)
9. start_targetSite: start location of the 3bp target locus
10. end_targetSite: end location of the 3bp target locus
11. read_type: Indicates if feature/gene is intronic, exonic, or intergenic.

NOTE: Depending on the information you decide to parse during your analysis, the information in this bam file may be different.  Throughout the explanation of the code logic below, we indicate the fields most important for analysis.  Please make changes to our scripts to ensure that the correct information is being processed if you decide to use another formatting for this bed file.


First we determine if a read has any indel in it by checking whether or not the CIGAR string indicates a perfect match.

In our experiment, we had reads spanning 76 nucleotides. Thus, if a CIGAR string read "76M", the read is not mutated, and is labeled as "WT"

*********************************
if ($4~/76M/){
	OFS="\t"
	# The read is a perfect match, label as WT
	print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"WT"
}
*********************************

NOTE: In our bed intersection files between reads and on-target sites, the 4th field contained the CIGAR string.  Change this if the CIGAR string is in a different field.  


If the read is NOT a perfect match, the CIGAR string is split into a vector where every other element in the vector is some sort of indel event. 

This is achieved using the following split command:

*********************************
split($4,a,/[[:digit:]].?M/,seps)
********************************* 

Example:
The CIGAR string 50M5D26M would be split into the following vector, called "a":
[BLANK,5D,BLANK]
The delimiters of this split are stored in another vector called "seps", and looks like this:
[50M,26M]

The indel vector is then processed in a for loop.  At the beginning of the for loop, the chromosomal locus of the read's start position is stored in the "strpos" variable.  When evaluating each indel site on the read, strpos is increased by the number of matching nucleotides preceeding the indel event:

********************************************************
strpos=strpos+substr(seps[i-1],1,(length(seps[i-1])-1))
********************************************************

The indel size is stored in the variable "msize":

**********************************
msize=substr(a[i],1,(length(a)-1))
**********************************


If the indel size is greater than the 3bp window where CRISPR mutations are expected to occur, the following logic is used to determine an indel event:

*******************************************
if (($9-strpos < msize) && ($9-strpos >=0))
*******************************************

Where field 9 is the start location of the expected CRISPR cut site.  (NOTE: change this field accordingly if this information is in a different field.)  


If the indel size is less than the 3bp window where CRISPR mutations are expected to occur, the following logic is used to determine an indel event:

**********************************************
if (($10-(strpos+msize)>=0) && (strpos-$9>=0))
**********************************************

Where field 10 is the end location of the expected CRISPR cut site. (NOTE: change this field accordingly if this information is in a different field.)

If either of the two logical statements are evaluated to be TRUE, "INDEL" is printed at the end of the line, the for loop is broken, and "INDEL" is printed at the end of the line.

However, if an indel event wasn't detected with both logical statements, the indel size is added to the "strpos" variable, and the rest of the indel sites from vector "a" are evaluated for the read.

If the loop manages to evaluate all the indel values in vector "a", then no indel events were found, and "WT" is printed at the end of the line.  

*************************************************************
if (i==length(a)){
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,"WT"
	}
*************************************************************


The script should be called like this:

awk -f indel_parse_whole.awk < On_target_reads.bed > On_target_reads_true.bed

The last field of "On_target_reads_true.bed" should have either "WT" or "INDEL". 


Using our "indel_freq_count.awk" script, users can create tables storing the number of legitimate CRISPR indels, wild-type sequences, and the fraction of indel events detected.  




















