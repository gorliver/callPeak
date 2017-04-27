#!/usr/bin/env perl
#===============================================================================
#
#         FILE: callPeak.pl
#        USAGE: ./callPeak.pl  
#  DESCRIPTION: call peak from DNase-seq data
#       AUTHOR: Hainan (), zhaohainancau@gmail.com
# ORGANIZATION: UW-Madison
#      VERSION: 1.0
#      CREATED: 04/25/2017 12:38:52 PM
#===============================================================================



use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);

use FindBin qw ($Bin);
use lib "$Bin/lib";

use smooth;
use threshold;


# get parameters
my $help = 0;
my $binSize=100;
my $stepSize=10;
my $outerBinSize=20000;
my $sigma=50;
my $randNum=50000;

my $bamFile;
my $chromFile;
my $uniqRegionFile;
my $bigWig="peak";
my $threads=4;

GetOptions(
    'help' => \$help,
    'b=s'   => \$bamFile,   # input bam file
	'g=s'   => \$chromFile, # chromosome length file
    'u=s'    => \$uniqRegionFile,   # unique regions file
    'o=s'   => \$bigWig,     # prefix of output 
    't=i'   => \$threads, # number of chromososme run parallelly
    's=i'    => \$stepSize,         # step size
    'l=i'    => \$binSize,          # bin size
	'bg=i'    => \$outerBinSize,          # background bin size
    'd=i'    => \$sigma,          # sigma
    'r=i'    => \$randNum,         # number of random regions for p-value calculation
) or pod2usage(2);

pod2usage(1) if $help;

if(!(defined $bamFile) || !(defined $chromFile) || !(defined $uniqRegionFile)){
    print STDERR "please set a bam file (-b), a chromosome file (-g) and a unique region file (-u).\n";
    pod2usage(1);
}

if(!(-e $bamFile)){
    print "\n***** $bamFile doesn't exist. *****\n\n";
    pod2usage(1);
}


if(!(-e $chromFile)){
    print "\n***** $chromFile doesn't exist. *****\n\n";
    pod2usage(1);
}

if(!(-e $uniqRegionFile)){
    print "\n***** $uniqRegionFile doesn't exist. *****\n\n";
    pod2usage(1);
}


# start calling peak
my $callPeak=smooth->new();

# generate pvalue table, get threthold

## make random regions

print "make random region\n";
print time,"\n";

$callPeak->setChromInfo($chromFile);
my @chr=map {$_->[0]} values %{$callPeak->getData('ChromInfo')};

print "chromosomes are got @chr\n";

$callPeak->makeRandPos($randNum,$binSize/2,$outerBinSize/2);

## make distribution
print "make distribution\n";
print time,"\n";
$callPeak->makeDistr($bamFile,$binSize/2,$outerBinSize/2,$uniqRegionFile);

## make pvala
print "generate pvalue table\n";
print time,"\n";
$callPeak->calPval();


## get pvalue
print "get pvalues\n";
print time,"\n";
my $pval=$callPeak->getData('Pvalue');

my $th;
foreach (sort {$a<=>$b} keys %{$pval}){
#    print "$_\t$$pval{$_}\n";
    $th=$_ if !(defined $th) && $$pval{$_}<0.01;
}

if(!defined $th){
    print "cannot determine the threthold!\n";
    exit;
}
else{
	print "threshold is $th\n";
}

$callPeak->gauss($binSize,$sigma);

my $fixed=$callPeak->getData('FixedArr');

my $thre=$th;
my @arr;
my @total;

my @test=randFeature($thre);
for(my $x=1;$x<=200;$x+=1){
    my $add=$test[$x] || 0;
    for(my $y=0;$y<@{$fixed};$y++){
        $arr[$y]+=$add*$$fixed[$y];
    }
    my @step=splice(@arr,0,$stepSize);
    my $sum;
    grep {$sum+=$_} @step;
    last if $sum==0;
    push @total,@step;
}


my $sd=sprintf("%.2f",stdev(\@total));

print "# sd is $sd\n";
print "#detect threshold\n";

# call peak, output bigwig and bed file
#
my $wigPm=Parallel::ForkManager->new($threads);
my @outFiles;


$wigPm->run_on_finish(
    sub {
        my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,$data_structure_reference ) =@_;
        if ( defined($data_structure_reference) ) {
            push @outFiles, $$data_structure_reference;
        }
        else{
            print STDERR "process $pid did not return a result, it probably failed for some reason. Run the script using single thread (-t 1) to see what happens\n\n";
        }
    }
);

WIG:
for(my $x=0;$x<@chr;$x++){
    my $pid=$wigPm-> start and next WIG;
    my $out=$bigWig . "$chr[$x]";
	print STDERR "analysing chr $chr[$x]\n";
    my $input;
    open $input,"samtools view -h $bamFile $chr[$x] |" or die;
    $callPeak->wigOut($input,$binSize,$outerBinSize,$stepSize,$uniqRegionFile,$out,$th,$sd);
    close $input;
    $wigPm->finish(0,\$out);
}
print STDERR "all job submitted\n";
$wigPm->wait_all_children;

print STDERR "all job finished,check files:\n@outFiles\n";


sub average{
    my $data = shift @_;
    if (not @$data) {
        die("Empty array\n");
    }
    my $total = 0;
    foreach (@$data) {
        $total += $_;
    }
    my $average = $total / @$data;
    return $average;
}
sub stdev{
    my $data = shift @_;
    if(@$data == 1){
        return 0;
    }
    my $average = average($data);
    my $sqtotal = 0;
    foreach(@$data) {
        $sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal / (@$data-1)) ** 0.5;
    return $std;
}

sub randFeature{
    my $thre=shift;
    my @arr;
    my @num;
    # get working array
    for(my $x=0;$x<=$thre;$x++){
        push @num,$x;
    }
    @num=(@num,reverse @num);
    # extract 20 or less element from working arr
    my $div=@num/20;
    my $re=@num%20;
    if($div<=1){
        @arr=@num;
    }
    else{
        my $ex=getEx($re,\@num);
        my $tt=int(@num/20);
        my $sw=0;
        for(my $x=0;$x<@num;$x+=$tt){
            push @arr,$num[$x] if !$$ex{$sw};
            $sw++;
        }
    }
    return @arr;
}

sub getEx{
    my $re=shift;
    my $numRef=shift;
    my $ex;
    my $sw=0;
    for(my $y=0;$y<100;$y++){
        my $tt=int(rand(@{$numRef})+0.5);

        if(!$$ex{$tt}){
            $$ex{$tt}++;
            $sw++;
        }

        last if $sw==$re;
    }
    return $ex;
}


__END__

=head1 NAME

callpeak.pl - call peak from DNase-seq data

=head1 SYNOPSIS

perl callPeak.pl -b <bam file> -g <chrom length file> -u <unique region file>

Call DNase-seq peaks using gausss kernel smooth. CallPeak.pl will normalize the read number by mappability and local background.

Options:

        -help             this help message.

          -b              input bam file, mandatory.

          -g              chromosome length file, mandatory.

                          This file contain two collomn: chromosome name and chromosome length:
                          chr1	300000000
                          chr2	280000000
                          ...

          -u              unique regions file, mandatory. This file need to be compressed by bgzip.
                          Please refer to http://www.htslib.org/doc/tabix.html for more information about bgzip.
                          This file contain three collomn: chromosome name, start position of the 
                          unique region and end position of the unique region:
                          chr1	1	100
                          chr1	200	300
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


=head1 DESCRIPTION

Call DNase-seq peaks using gausss kernel smooth. CallPeak.pl will normalize the read number by mappability and local background.

=cut


