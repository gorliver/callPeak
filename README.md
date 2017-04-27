# callPeak
call peaks from DNase-seq data


# requirement

    perl
    perl module:
    Parallel::ForkManager
    Bio::DB::HTS::Tabix

    python
    bigWigpy (https://github.com/dpryan79/pyBigWig)

    samtools
    bedtools

# installation

## 1. Install perl modules:
Using cpanm is probably the easiest why to install perl modules.

    ./cpanm Parallel::ForkManager
    ./cpanm Bio::DB::HTS::Tabix

    

## 2. Install pyBigWig
Please refer to https://github.com/dpryan79/pyBigWig.

## 3. install samtools
Please refer to http://www.htslib.org/ for instruction
    
## 3. install callPeak
    git clone https://github.com/gorliver/callPeak
    cd callPeak
    chmod +x callPeak.pl

Or simplely copy the whole directory to you local mathine.
    
Then add the path of callPeak.pl to PATH.

# usage

## 1. generate uniquely mapped regions.
To reduce the influence of repeat sequences, callPeak.pl will normalize the read number by mappability of each bins.

### 1. simulate artifical reads from the target genome. 
The read length of the simulated reads should be the same as the real data. For example, simulating single end reads with 20bp length in Arabidopsis genome:
  
    wgsim -1 20 -h -N 2000000 Tair10.fa out.read1.fq /dev/null
    
### 2. Mapped the simulated reads to target genome.
Mapping the simulated reads to target genome using your favorate aligner. In the following command, only reads with MQ>=20 are retained.
    
    bwa aln Tair10.fa out.read1.fq | bwa samse Tair10.fa - out.read1.fq | samtools view -q 20 -bS - | samtools sort - -o simu_sort.bam
    
### 3. extract uniquely mapped regions
    bamToBed -i simu_sort.bam | mergeBed -d 1 -i - > uniqRegionAra20bp.bed

### 4. compressed the file using bgzip

    bgzip uniqRegionAra20bp.bed
    
For more information about bgzip, please check http://www.htslib.org/doc/tabix.html
    
## 2. run callPeak.pl    
    callPeak.pl Tair10DNase-seq_sort.bam chromAra.len uniqRegionAra20bp.bed.gz araPeak 5

This commond will call peaks from bam file. Note that the input bam file need to be sorted.
    
    
    
Here is the description of the options:

    perl callPeak.pl -b <bam file> -g <chrom length file> -u <unique region file>

    Call DNase-seq peaks using gausss kernel smooth. CallPeak.pl will
    normalize the read number by mappability and local background.

    Options:

            -help             this help message.

              -b              input bam file, mandatory.

              -g              chromosome length file, mandatory.

                              This file contain two collomn: chromosome name and chromosome length:
                              chr1  300000000
                              chr2  280000000
                              ...

              -u              unique regions file, mandatory. This file need to be compressed by bgzip.
                              Please refer to http://www.htslib.org/doc/tabix.html for more information about bgzip.
                              This file contain three collomn: chromosome name, start position of the
                              unique region and end position of the unique region:
                              chr1  1       100
                              chr1  200     300
                              ...

              -o              prefix of output, default is "peak".

              -t              number of chromososme run parallelly.

                              for example, if there are 5 chromsomes in the bam file and you have at
                              least 5 threads in your machine, then set "-t 5" to run the 5 chromosomes
                              at same time. "-t 1 " means callpeak.pl will analyze each chromsome one by one.

              -s              step size, default is 10. Smaller value will increase running time but give
                              slightly better smoothness.

              -l              bin size. default is 100. Smaller value will give sharper peak, larger value
                              will give better smoothness.

              -bg             background bin size. default is 20000. Used to estimate the local backgroud
                              around each bin. Don't set too small, otherwise the broad peak well be missed.

              -d              sigma for kernel smooth. default is 50.

              -r              number of random regions for p -value calculation. default is 50000.

