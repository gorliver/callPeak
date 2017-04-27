# callPeak
call peaks from DNase-seq data


# requirement

    perl
    perl module:
    Parallel::ForkManager
    Bio::DB::HTS::Tabix

    python
    bigWigpy (https://github.com/dpryan79/pyBigWig)

# installation

## 1. Install perl modules:
    ./cpanm Parallel::ForkManager
    ./cpanm Bio::DB::HTS::Tabix

    cpanm is provided in the directory

## 2. Install pyBigWig
    Please refer to https://github.com/dpryan79/pyBigWig.

## 3. install callPeak
    git clone https://github.com/gorliver/callPeak
    cd callPeak
    chmod +x callPeak.pl

    Then add the path of callPeak.pl to PATH.

# usage

    perl callPeak.pl -gs <genome file> -t <target sequence file>

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

