# Marine Biology 740 Short Project 2

Group: Oscar Ramfelt, 

Link to Github repository [here](https://github.com/kolaban/mbio740-short-project2)

## Download

Shotgun sequencing data was downloaded from SRX6892810 using the following command:

```Bash
fastq-dump --split-files SRX6892810
```

Which downloaded the forward (SRR10167988_1.fastq) and reverse (SRR10167988_2.fastq) reads. We also downloaded the reference genome for *Staphylococcus aureus* using the ftp link to its assembly found [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/425/GCA_000013425.1_ASM1342v1/GCA_000013425.1_ASM1342v1_genomic.fna.gz) (GCA_000013425.1_ASM1342v1_genomic.fna.gz).

### Raw Coverage

We found the coverage of the reference genome is the shotgun sequencing data by using bowtie2 (version = 2.3.5.1). The two following command was used when running bowtie2:

```Bash
bowtie2-build GCA_000013425.1_ASM1342v1_genomic.fna.gz staphylococcus
```

```Bash
bowtie2 --threads 6 -x reference/staphylococcus \
        --no-unal \
        -1 SRR10167988_1.fastq \
        -2 SRR10167988_2.fastq \
        -S staphylococcus_bowtie.sam
```

Following this we then used samtools to find the average coverage depth.

```Bash
samtools view -F 4 -bS staphylococcus_bowtie.sam > staphylococcus_bowtie-RAW.bam
samtools sort staphylococcus_bowtie-RAW.bam -o staphylococcus_bowtie.bam
samtools index staphylococcus_bowtie.bam
```

Calculation of coverage was done using samtools and based off the helpful advice of [this](https://sarahpenir.github.io/bioinformatics/awk/calculating-mapping-stats-from-a-bam-file-using-samtools-and-awk/) website.

```Bash
samtools depth -a staphylococcus_bowtie.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'
```

This gave us a coverage of **93.4082%** across the genome.

## Fastqc

### Output

The raw reads were then quality checked using fastqc through the following command:

```Bash
fastqc SRR10167988_1.fastq
fastqc SRR10167988_2.fastq
```

Producing two html files with the reads associated statistics, the html files can be found in the read_quality directory. Based on these two html reports it seems that overall the quality of the data is good with some minor issues. The graphs were roughly i dentical between the forward and reverse reads so only the graphs for forward reads.

#### Per Base Sequence Quality

![Fig 1](static/before_qc/per_base_seq_quality.png)

#### Per Base Sequence Content

![Fig 2](static/before_qc/seq_content_issue.png)

#### Sequence Duplication Levels

![Fig 3](static/before_qc/seq_duplication_issue.png)

#### Adapter Content

![Fig 4](static/before_qc/adapter_issue.png)

### Cleaning Strategy

The overall cleaning strategy will be to remove the first few bases since they are both lower quality compared to the others and show a strange pattern in the per base sequence content graph. It also seems that there is a signficant amount of adapter content on the tail end of the reads showing that any cleaning strategy should include adapter trimming to remove those issues.

## Trimmomatic

To clean up the sequences we used Trimmomatic (version = 0.39), we decided to remove the first few bases and the last two to reduce the issues that Fastqc showed us as well as add an adapter trimming step to remove the adapter content that Fastqc showed us was an issue. Based on the Fastqc graph it seems that the adapters are Nextera adapters so we used the pair end adapters that come default with trimmomatic i.e. NexteraPE-PE.fa.

```Bash
trimmomatic PE SRR10167988_1.fastq SRR10167988_2.fastq \
            forward_good forward_fail reverse_good reverse_fail \
            ILLUMINACLIP:/opt/conda/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
            LEADING:3 TRAILING:2
```

### Post QC Results

Despite the previous efforts to do qc with trimmomatic it seems that the per base sequence content remains an issue. However, after some reseach online about the nextera library prep method it seems that this is likely due to the library prep itself and one can choose whether or not to remove it. More information found [here](http://seqanswers.com/forums/showthread.php?t=45024).

There were some warnings in the data that did not show up previously but they did not seem too severe and were likely due to the fact that we had trimmed off sequences. The three graphs where we had warnings/fails are shown below.

#### Per Base Sequence Content

![Fig 6](static/after_qc/seq_content_issue.png)

#### Sequence Length Distribution

![Fig 5](static/after_qc/seq_length_issue.png)

#### Sequence Duplication Levels

![Fig 6](static/after_qc/seq_duplication_issue.png)

In addition to running Fastqc we also ran bowtie2 again to check for any change in coverage across the *Staphylococcus aureus* genome. Our command is identical to what we ran earlier except that we are now using the cleaned up fastq files.

```Bash
bowtie2 --threads 6 -x ../reference/staphylococcus \
        --no-unal \
        -1 ../qc_reads/forward_good \
        -2 ../qc_reads/reverse_good \
        -S staphylococcus_bowtie_post.sam
```

We then used the same method outlined previously to calculate average coverage.

```Bash
samtools view -F 4 -bS staphylococcus_bowtie_post.sam > staphylococcus_bowtie_post-RAW.bam
samtools sort staphylococcus_bowtie_post-RAW.bam -o staphylococcus_bowtie_post.bam
samtools index staphylococcus_bowtie_post.bam
```

Coverage calculations also followed the same method.

```Bash
samtools depth -a staphylococcus_bowtie_post.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'
```

Giving use a coverage of **93.3981%** which is almost no change from our estimation prior to cleaning, which probably means that none of the reads that we removed were core to the genome that we want to assemble.

## Spades

For the assembly we used Spades (version = ) with the following command being run. 

```Bash

```
