# RNA-seq-2

RNA-seq analysis in R
Obtaining and aligning RNA-seq reads from public repositories
Stephane Ballereau, Mark Dunning, Oscar Rueda, Ashley Sawle
Last modified: 26 Apr 2017
Set up a database of SRA runs
Raw reads from NGS experiments tend to be distributed through the Short Read Archive (SRA). The SRAdb Bioconductor package can be used to query and download files that are hosted in SRA. More information can be found in the package vignette.

vignette("SRAdb")
Firstly, we need to download a database file. This is a large file and could take some time (25.82 )Gb.

library(SRAdb)
sqlfile <-'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
Obtain information for a particular experiment
We can now query what information is available for a particular experiment; in this case SRP045534. This should list the samples that are available and their respective identifiers.

sraInf <- getSRAinfo("SRP045534",sra_con, sraType="sra")
sraInf
Download the set of sra files
We can now directly download each sra file. The sra file is SRA’s own archive format, but we can extract the raw reads in the more common .fastq format in the next step.

Here, sapply is a convenient way of repeating the same operation for a vector of arguments. In this case we want to run the getSRAfile function with different file names.

sapply(sraInf$run, function(x) try(getSRAfile(x,sra_con, fileType="sra"),silent=TRUE))
Extracting fastq files
Using the sra-toolkit command-line utility from NCBI we can generate the fastq files from these archive files. We can do this within a Terminal (i.e. not within RStudio) with the following, making sure your working directory contains the .sra files.

for sra in *.sra
do
fastq-dump $sra
done
After each fastq file has been extracted, you should see a message to report have many reads (spots) are contained in the file

Quality assessment of reads
The fastqc is recommened for a preliminary assessment of the read quality. However, caution should be exercised when interpreting the results as the reports are not specifically-tailored for RNA-seq. Some sections are know to flag-up warning or error messages for perfectly fine RNA-seq experiments.

We can run this at the command line:-

for fq in *.fastq
do
fastqc $fq
done
Downloading the reference genome
To align to the mm10 genome, we will first download the reference genome from UCSC. We have to download each individual chromosome separately, and then join together into a single file.

wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
gunzip chromFa.tar.gz
tar xvf chromFa.tar
cat *.fa > mm10.fa
rm chr*.fa
rm chromFa.tar.gz
Alignment using bowtie
Firstly, we need to build an index file from the reference genome that we have downloaded:-

bowtie2-build mm10.fa mm10
In practice, we would probably run the alignment of each sample in parallel using the high-performance cluster. However, for illustration purposes, we give the script that will align each sample individually.

bowtie2 -x mm10 -U SRR1552444.fastq -S SRR1552444.sam
samtools view -bS SRR1552444.sam > SRR1552444.bam
samtools sort SRR1552444.bam -o SRR1552444.sorted.bam
samtools index SRR1552444.sorted.bam

bowtie2 -x mm10 -U SRR1552445.fastq -S SRR1552445.sam
samtools view -bS SRR1552445.sam > SRR1552445.bam
samtools sort SRR1552445.bam -o SRR1552445.sorted.bam
samtools index SRR1552445.sorted.bam

bowtie2 -x mm10 -U SRR1552446.fastq -S SRR1552446.sam
samtools view -bS SRR1552446.sam > SRR1552446.bam
samtools sort SRR1552446.bam -o SRR1552446.sorted.bam
samtools index SRR1552446.sorted.bam

bowtie2 -x mm10 -U SRR1552447.fastq -S SRR1552447.sam
samtools view -bS SRR1552447.sam > SRR1552447.bam
samtools sort SRR1552447.bam -o SRR1552447.sorted.bam
samtools index SRR1552447.sorted.bam

bowtie2 -x mm10 -U SRR1552448.fastq -S SRR1552448.sam
samtools view -bS SRR1552448.sam > SRR1552448.bam
samtools sort SRR1552448.bam -o SRR1552448.sorted.bam
samtools index SRR1552448.sorted.bam

bowtie2 -x mm10 -U SRR1552449.fastq -S SRR1552449.sam
samtools view -bS SRR1552449.sam > SRR1552449.bam
samtools sort SRR1552449.bam -o SRR1552449.sorted.bam
samtools index SRR1552449.sorted.bam

bowtie2 -x mm10 -U SRR1552450.fastq -S SRR1552450.sam
samtools view -bS SRR1552450.sam > SRR1552450.bam
samtools sort SRR1552450.bam -o SRR1552450.sorted.bam
samtools index SRR1552450.sorted.bam

bowtie2 -x mm10 -U SRR1552451.fastq -S SRR1552451.sam
samtools view -bS SRR1552451.sam > SRR1552451.bam
samtools sort SRR1552451.bam -o SRR1552451.sorted.bam
samtools index SRR1552451.sorted.bam

bowtie2 -x mm10 -U SRR1552452.fastq -S SRR1552452.sam
samtools view -bS SRR1552452.sam > SRR1552452.bam
samtools sort SRR1552452.bam -o SRR1552452.sorted.bam
samtools index SRR1552452.sorted.bam

bowtie2 -x mm10 -U SRR1552453.fastq -S SRR1552453.sam
samtools view -bS SRR1552453.sam > SRR1552453.bam
samtools sort SRR1552453.bam -o SRR1552453.sorted.bam
samtools index SRR1552453.sorted.bam

bowtie2 -x mm10 -U SRR1552454.fastq -S SRR1552454.sam
samtools view -bS SRR1552454.sam SRR1552454.bam
samtools sort SRR1552454.bam -o SRR1552454.sorted.bam
samtools index SRR1552454.sorted.bam

bowtie2 -x mm10 -U SRR1552455.fastq -S SRR1552455.sam
samtools view -bS SRR1552455.sam > SRR1552455.bam
samtools sort SRR1552455.bam -o SRR1552455.sorted.bam
samtools index SRR1552455.sorted.bam
Renaming to be consistent with GEO
The files we have just created are named according to their SRA identifier. However, these names are not very useful for analysis. The Gene Expression Omnibus (GEO) entry for the dataset has the mapping information between SRA and sample identifers.

library(GEOquery)
tmp <- getGEO("GSE60450")
gseInf <- pData(tmp[[1]])
gseInf
 
 
title
<fctr>
geo_accession
<fctr>
status
<fctr>
submission_date
<fctr>
GSM1480291	Luminal virgin #1	GSM1480291	Public on Jan 19 2015	Aug 15 2014	
GSM1480292	Luminal virgin #2	GSM1480292	Public on Jan 19 2015	Aug 15 2014	
GSM1480293	Luminal 18.5 dP #1	GSM1480293	Public on Jan 19 2015	Aug 15 2014	
GSM1480294	Luminal 18.5 dP #2	GSM1480294	Public on Jan 19 2015	Aug 15 2014	
GSM1480295	Luminal 2 dL #1	GSM1480295	Public on Jan 19 2015	Aug 15 2014	
GSM1480296	Luminal 2 dL #2	GSM1480296	Public on Jan 19 2015	Aug 15 2014	
GSM1480297	Basal virgin #1	GSM1480297	Public on Jan 19 2015	Aug 15 2014	
GSM1480298	Basal virgin #2	GSM1480298	Public on Jan 19 2015	Aug 15 2014	
GSM1480299	Basal 18.5 dP #1	GSM1480299	Public on Jan 19 2015	Aug 15 2014	
GSM1480300	Basal 18.5 dP #2	GSM1480300	Public on Jan 19 2015	Aug 15 2014	
1-10 of 12 rows | 1-5 of 52 columns
We obtain a new name for each bam file by joining the metadata from SRA and GEO.

library(dplyr)
sraInf <- mutate(sraInf, bam=paste0(run, ".sorted.bam"))
gseInf <- mutate(gseInf, experiment = basename(as.character(supplementary_file_2)),
                 newbam = gsub("Sample name: ","", description),
                 newbam = gsub("-",".",newbam,fixed=TRUE),
                 newbam = paste0(newbam, ".bam"))
gseInf
title
<fctr>
geo_accession
<fctr>
status
<fctr>
submission_date
<fctr>
Luminal virgin #1	GSM1480291	Public on Jan 19 2015	Aug 15 2014	
Luminal virgin #2	GSM1480292	Public on Jan 19 2015	Aug 15 2014	
Luminal 18.5 dP #1	GSM1480293	Public on Jan 19 2015	Aug 15 2014	
Luminal 18.5 dP #2	GSM1480294	Public on Jan 19 2015	Aug 15 2014	
Luminal 2 dL #1	GSM1480295	Public on Jan 19 2015	Aug 15 2014	
Luminal 2 dL #2	GSM1480296	Public on Jan 19 2015	Aug 15 2014	
Basal virgin #1	GSM1480297	Public on Jan 19 2015	Aug 15 2014	
Basal virgin #2	GSM1480298	Public on Jan 19 2015	Aug 15 2014	
Basal 18.5 dP #1	GSM1480299	Public on Jan 19 2015	Aug 15 2014	
Basal 18.5 dP #2	GSM1480300	Public on Jan 19 2015	Aug 15 2014	
1-10 of 12 rows | 1-4 of 54 columns
combinedInf <- left_join(gseInf, sraInf, by="experiment")
combinedInf %>% select(description,description.1,experiment,bam,newbam)
description
<fctr>
description.1
<fctr>
experiment
<chr>
Sample name: MCL1-LA	MCL1-LA_BC2CTUACXX_GATCAG_L001_R1	SRX681985	
Sample name: MCL1-LB	MCL1-LB_BC2CTUACXX_TGACCA_L001_R1	SRX681986	
Sample name: MCL1-LC	MCL1-LC_BC2CTUACXX_GCCAAT_L001_R1	SRX681987	
Sample name: MCL1-LD	MCL1-LD_BC2CTUACXX_GGCTAC_L001_R1	SRX681988	
Sample name: MCL1-LE	MCL1-LE_BC2CTUACXX_TAGCTT_L001_R1	SRX681989	
Sample name: MCL1-LF	MCL1-LF_BC2CTUACXX_CTTGTA_L001_R1	SRX681990	
Sample name: MCL1-DG	MCL1-DG_BC2CTUACXX_ACTTGA_L002_R1	SRX681991	
Sample name: MCL1-DH	MCL1-DH_BC2CTUACXX_CAGATC_L002_R1	SRX681992	
Sample name: MCL1-DI	MCL1-DI_BC2CTUACXX_ACAGTG_L002_R1	SRX681993	
Sample name: MCL1-DJ	MCL1-DJ_BC2CTUACXX_CGATGT_L002_R1	SRX681994	
1-10 of 12 rows | 1-3 of 5 columns
combinedInf
title
<fctr>
geo_accession
<fctr>
status
<fctr>
submission_date
<fctr>
Luminal virgin #1	GSM1480291	Public on Jan 19 2015	Aug 15 2014	
Luminal virgin #2	GSM1480292	Public on Jan 19 2015	Aug 15 2014	
Luminal 18.5 dP #1	GSM1480293	Public on Jan 19 2015	Aug 15 2014	
Luminal 18.5 dP #2	GSM1480294	Public on Jan 19 2015	Aug 15 2014	
Luminal 2 dL #1	GSM1480295	Public on Jan 19 2015	Aug 15 2014	
Luminal 2 dL #2	GSM1480296	Public on Jan 19 2015	Aug 15 2014	
Basal virgin #1	GSM1480297	Public on Jan 19 2015	Aug 15 2014	
Basal virgin #2	GSM1480298	Public on Jan 19 2015	Aug 15 2014	
Basal 18.5 dP #1	GSM1480299	Public on Jan 19 2015	Aug 15 2014	
Basal 18.5 dP #2	GSM1480300	Public on Jan 19 2015	Aug 15 2014	
1-10 of 12 rows | 1-4 of 61 columns
The base R function file.symblink can be used to create symbolic links from one file to another; thus retaining the original file name and avoid creating a complete copy of each file. Such links are often used in NGS data when we don’t want to create copies of files that are potentially rather large. With this approach, when we want to access MCL1.LA.bam (for example), the file system will know to actually access SRR1552444.sorted.bam.

for(i in seq_along(combinedInf$bam)){
  
  file.symlink(combinedInf$bam[i], combinedInf$newbam[i])
  file.symlink(paste0(combinedInf$bam[i],".bai"), paste0(combinedInf$newbam[i],".bai"))
  
}
list.files()
Alignment using Rsubread
Alignment could also be performed using Rsubread as we did in the first practical session. The only difference is to use the entire mm10 genome and point to the newly-downloaded fastq files.

library(Rsubread)
buildindex("mm10",reference="mm10.fa")
fastqfiles <- list.files(pattern=".fastq")
align("mm10",readfile1=fastqfiles)
