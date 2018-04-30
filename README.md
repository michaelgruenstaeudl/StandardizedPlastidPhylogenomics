Scripts for Standardized Plastid Phylogenomics
==============================================
*Version: 2018.04.30, Author: Michael Gruenstaeudl*

Software scripts to automate and standardize the various bioinformatic processes in plastid phylogenomics


GETTING STARTED
---------------
##### File prerequisites
* Raw reads from Illumina platform (R1 and R2, in separate files)
* Control and sample number of your Illumina run (i.e., the last two fields of a FASTQ header)
  * Note: For more info the format of a FASTQ header, see [here](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm)
* One or more plastid genomes that can be used as reference genomes (in a single FASTA-file)

##### Software prerequisites
The following software must available on your system:
* [Git](https://git-scm.com/) (for installation)
* [bioawk](https://github.com/lh3/bioawk) (for script 1)
* [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html) (for script 2)
* [Python 2.7](https://www.python.org/download/releases/2.7/) (for scripts 3,7)
* [IOGA](https://github.com/holmrenser/IOGA) (for script 3)
* [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (for script 4)
* [bowtie2](http://bowtie-bio.sf.net/bowtie2) (for scripts 5,6)
* [samtools](https://github.com/samtools/samtools) (for scripts 5,6)
* [bedtools](https://github.com/arq5x/bedtools2) (for scripts 5,6)
* [Perl 5](https://www.perl.org/) (for script 7)

INSTALLATION
------------
##### Cloning this GitHub repository
```
cd /path_to_git/
git clone https://github.com/michaelgruenstaeudl/StandardizedPlastidGenomeAssemblyAnnotation.git
```

USAGE
-----

### 1. Specifying file names and locations
##### Specifying the working directory and the local Git repository
```
MYWD=/path_to_infiles/
MYGIT=/path_to_git/StandardizedPlastidGenomeAssemblyAnnotation/
```
##### Specifying sample name, reference name and FASTQ header
```
MYSAMPLE="Sample_name_here"  # Example: MYSAMPLE="NY684"
MYREF="Reference_name_here"  # Example: MYREF="MG720559"
FASTQ_HEADR="controlNum_sampleNum_here"  # Example: FASTQ_HEADR="N:0:6"
```
##### Specifying and initializing LOG-file
Adding the project name as the first line of the log-file
```
LOG=$MYWD/${MYSAMPLE}.log
echo "$MYSAMPLE" > $LOG
```

### 2. Quality filtering of reads
##### SCRIPT 1 - Generating ordered intersection of reads
Bash script to generate an ordered intersection of paired-end Illumina reads and to quantify the loss of reads during this process
```
## DECLARING VARIABLES
INF1=$MYWD/${MYSAMPLE}_R1.fastq
INF2=$MYWD/${MYSAMPLE}_R2.fastq

## RUNNING SCRIPT
bash $MYGIT/Script1.sh $INF1 $INF2 $LOG $FASTQ_HEADR
```
##### SCRIPT 2 - Conducting quality filtering of reads
Bash script to conduct quality filtering of Illumina reads and to quantify the loss of reads during quality filtering.
```
## DECLARING VARIABLES
INF1=$MYWD/${MYSAMPLE}_R1.intersect.fastq
INF2=$MYWD/${MYSAMPLE}_R2.intersect.fastq

## RUNNING SCRIPT
bash $MYGIT/Script2.sh $INF1 $INF2 $LOG
```

### 3. Plastid genome assembly
##### SCRIPT 3 - Conducting the genome assembly process
Bash script to automate the genome assembly process and to quantify assembly success.
```
## DECLARING VARIABLES
INF1=$MYWD/${MYSAMPLE}_R1.Q30filt.fastq
INF2=$MYWD/${MYSAMPLE}_R2.Q30filt.fastq
IOGA=/path_to_IOGA/IOGA.py
REF_DB=/path_to_reference_genomes/combined_references.fas

## RUNNING SCRIPT
bash $MYGIT/Script3.sh $INF1 $INF2 $LOG $IOGA $REF_DB $MYSAMPLE
```
##### Manual step - Stitching of contigs to complete genome
Recommendation: Use [Geneious](https://www.geneious.com/) for this step


### 4. Evaluation of assembly
##### SCRIPT 4 - Confirming the IR boundaries
Bash script to automate the confirmation of IR boundaries via self-blasting.
```
## DECLARING VARIABLES
FINAL_ASMBLY=$MYWD/${MYREF}.fasta

## RUNNING SCRIPT
bash $MYGIT/Script4.sh $FINAL_ASMBLY $LOG
```
##### SCRIPT 5 - Extracting reads that map to final assembly, Generating assembly statistics - part 1
Bash script to automatically map the quality-filtered reads against the final assembly and to extract the mapped paired reads from the quality-filtered read files. This script also compiles a series of mapping statistics to describe the mapping process.
```
## DECLARING VARIABLES
INF1=$MYWD/${MYSAMPLE}_R1.Q30filt.fastq
INF2=$MYWD/${MYSAMPLE}_R2.Q30filt.fastq
FINAL_ASMBLY=$MYWD/${MYREF}.fasta

## RUNNING SCRIPT
bash $MYGIT/Script5.sh $INF1 $INF2 $FINAL_ASMBLY $LOG $MYSAMPLE
```
##### SCRIPT 6 - Generating assembly statistics - part 2
Bash script to automatically generate assembly quality statistics (i.e., the mean read length of mapped reads and the number of nucleotides with coverage depth equal or greater than 20, 50 and 100)
```
## DECLARING VARIABLES
INF_MAPPED_R1=$MYWD/${MYSAMPLE}.MappedAgainst.${MYREF}_R1.fastq
INF_MAPPED_R2=$MYWD/${MYSAMPLE}.MappedAgainst.${MYREF}_R2.fastq
INF_SAM=$MYWD/${MYSAMPLE}.MappedAgainst.${MYREF}.sam

## RUNNING SCRIPT
bash $MYGIT/Script6.sh $INF_MAPPED_R1 $INF_MAPPED_R2 $INF_SAM $LOG
```

### 5. Plastid genome annotation

##### Manual step - Plastid genome annotation via DOGMA and cpGAVAS
* Annotation server [DOGMA](http://dogma.ccbb.utexas.edu/)
    * Note: Save the annotations of DOGMA as "Plain Text Summary"
* Annotation server [cpGAVAS](http://www.herbalgenomics.org/cpgavas/)

##### SCRIPT 7 - Convert DOGMA plain text summary to GFF
Bash script to convert the plain text summary of the annotations generated by DOGMA to a GFF file.
```
## DECLARING VARIABLES
INF=$MYWD/${MYSAMPLE}_DOGMA_Annotations_PlainTextSummary.txt
FINAL_ASMBLY=$MYWD/${MYREF}.fasta

## RUNNING SCRIPT
bash $MYGIT/Script7.sh $INF $FINAL_ASMBLY
```

##### SCRIPT 8 - Combine annotations into single set and generate union set of the annotations
Bash script to combine the annotations produced by the annotation servers DOGMA and cpGAVAS and to generate a union set of the combined annotations
```
## DECLARING VARIABLES
INF_CPGAVAS=$MYWD/${MYSAMPLE}_cpGAVAS_Annotations.gff
INF_DOGMA=$MYWD/${MYSAMPLE}_DOGMA_Annotations_PlainTextSummary.gff
FINAL_ASMBLY=$MYWD/${MYREF}.fasta
IRB="start_of_IRb,end_of_IRb"  # Example: IRB="89435,114903"
IRA="start_of_IRa,end_of_IRa"  # Example: IRA="134019,159487"

## RUNNING SCRIPT
bash $MYGIT/Script8.sh $INF_CPGAVAS $INF_DOGMA $FINAL_ASMBLY $IRB $IRA
```

##### SCRIPT 9 - Foo bar baz
Bash script to Foo bar baz
```
## Foo bar baz
Foo bar baz

## RUNNING SCRIPT
bash $MYGIT/Script9.sh Foo bar baz
```

##### SCRIPT 10 - Foo bar baz
Bash script to Foo bar baz
```
## Foo bar baz
Foo bar baz

## RUNNING SCRIPT
bash $MYGIT/Script10.sh Foo bar baz
```

##### SCRIPT 11 - Foo bar baz
Bash script to Foo bar baz
```
## Foo bar baz
Foo bar baz

## RUNNING SCRIPT
bash $MYGIT/Script11.sh Foo bar baz
```

TEST DATA
---------
Test data is available at: https://zenodo.org/record/1213269

CITATION
--------
Gruenstaeudl M., Gerschler N., Borsch T. (2018) Sharing bioinformatic workflows for plastid genomics - An example from the sequencing of complete plastid genomes of *Cabomba* (Cabombaceae). In Review.


CHANGELOG
---------
* 2018.04.30 - Update
   * Addition of scripts 9, 10, 11
* 2018.04.05 - Update
   * Addition of script 7, update of script 8
* 2018.04.04 - Update
   * Update of README
* 2018.04.03 - Update
   * Update of scripts
* 2018.02.28 - Original Upload
   * First upload of scripts
