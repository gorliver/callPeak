package smooth;

use strict;
use warnings;
use parent qw(threshold);


sub new{
	my $class=shift;
	my $self={};
	return bless $self,$class;
}

sub gauss{
	my $self=shift;
	my $binSize=shift;
	my $sigma=shift;
	my $st=-1*$binSize/2;
	my $ed=1*$binSize/2;
	my $fixedArr;

	for(my $x=$st;$x<=$ed;$x++){
		my $gsout=$self->kernel($x,$sigma);
		push @{$fixedArr},$gsout/$sigma; # add /$sigma
	}

	$self->{'FixedArr'}=$fixedArr;
}

sub getFixedArr{
	my $self=shift;
	return $self->{'FixedArr'};
}

sub kernel{
	my $self=shift;
    my $par=shift;
	my $sigma=shift;
    my $ratiosq=($par/$sigma)**2;
    my $cont=1/(sqrt(2*3.1415926));
    my $gs=$cont*exp(-0.5*($ratiosq));
    $gs=sprintf("%.5f",$gs);
    return $gs;
}

sub wigOut{
	my $self=shift;
	my $bamFile=shift;
	my $innerBinSize=shift;
	my $outerBinSize=shift;
	my $stepSize=shift;
	my $uniqRegionFile=shift;
	my $outWig=shift;
	my $cutoff=shift || 9;
	my $bound=shift || 0.2;
	my $arrInnerRef;
	my $arrOuterRef;
	my $innerBinStart=1;
	my $innerBinEnd=$innerBinStart+$innerBinSize;
	my $outerBinStart=1;
	my $outerBinEnd=$outerBinStart+$outerBinSize;
	my $Nchr="NA";
	my $densityRef;
	my $dbg=0;
	my $gap=0;
	
	my $FixedArr=$self->getFixedArr();

	my $trackName="$bamFile.$innerBinSize\.$stepSize";
	my $description=$trackName;
	my $bw=Inline::Python::Object->new('__main__','bwObj');
#	print "## create a bigwig file\n";
#	$bw->openBw($outWig,"w");
#	print "## add header\n";
#	$bw->addHead($chr,150982314);
	my $bwStart=1;
	my @bwEntries;

	my $mapped=$self->getmapReg($uniqRegionFile);
	my @mappedSt;
	my @mappedEd;

	my $peakTh=0;
	my $peakOn=0;
	my $peakStart=0;
	my $peakEnd=0;
	my $peakScore=0;

	my $dens=0;
	my $dens1=0;
	my $dens2=0;
	my $dens3=0;

	my $dis=0;
	
	my $tt=time;
	my $aa=1;
	
	my $genomeSize=0;
	my %chrLen;

	my $pval=$self->getData('Pvalue');

	while(<$bamFile>){
		if(/^\@/){
			while(<$bamFile>){
				last if !/^\@/;
				next if /sca/;
				next if /Mt/;
				next if /Pt/;
				if(/SN:(\S+)\s+LN:(\d+)/){$genomeSize+=$2;$chrLen{$1}=$2}
			}
		}
		die("no header in bam/sam file\n") if $genomeSize==0;
		
		my ($chr,$flag,$pos,$qual,$tag)=(split)[2,1,3,4,5];
		next if $qual<20;
		next if !$chrLen{$chr};
#		next if $chr eq "Mt" || $chr eq "Pt";
#		next if $chr==10 || $chr==9 || $chr==8;
#		next if ($chr==6 && $pos<150000000);
#		next if $chr!~/\d/; # need edit
		my $cutPos=$pos;
		if($flag==16){ # get cutting position of reverse reads
			grep {$cutPos+=$_} ($tag=~/(\d+)[MD]/g);
		}

#		print "$chr\t$pos\t",$pos+1,"\tTAG\n";
		if($chr ne $Nchr){#reach new chr or start from the first chr;
    		print STDERR "new chr ## add header $chr, $Nchr\n";
			if(@bwEntries>0){ # output the last calculation
					print "we got 3 $bwStart\n";
					$bw->addEnt($Nchr,$bwStart,\@bwEntries,1,1);
					undef @bwEntries;
					$bw->close();
			}
#			$bw->close() if -e "$outWig\_$Nchr";

			$Nchr=$chr;
			my $outputFile="$outWig\_$Nchr\.bw";
			my $outputBed="$outWig\_$Nchr\.bed";
			open BED,"+>$outputBed" or die ("cannot create bed file\n");
			$bw->openBw($outputFile,"w");
			$bw->addHead($Nchr,$chrLen{$Nchr});

			print STDERR "Start Chr $Nchr\n";

			$innerBinStart=1;
			$innerBinEnd=$innerBinStart+$innerBinSize;
			$outerBinStart=1;
			$outerBinEnd=$innerBinEnd+($outerBinSize/2-$innerBinSize/2);#$outerBinStart+$outerBinSize;

			undef $densityRef;
			undef $arrInnerRef;
			undef $arrOuterRef;

			$bwStart=$innerBinStart;

#			if(@bwEntries>0){ # output the last calculation
#				print "we got 3 $bwStart\n";
#				$bw->addEnt($Nchr,$bwStart,\@bwEntries,1,1);
#				undef @bwEntries;
#			}
#

			$peakTh=0;
			$peakOn=0;
			$peakEnd=0;
			$dens=0;
			$dens1=0;
			$dens2=0;
			$dens3=0;
			@mappedSt=sort {$a<=>$b} keys %{$$mapped{$Nchr}};
			@mappedEd=sort {$a<=>$b} values %{$$mapped{$Nchr}};

			redo;
		}

		($tt,$aa)=timer($pos, $aa, $tt, $genomeSize);

	
		push @{$arrInnerRef},$cutPos; # record the cuttings
		push @{$arrOuterRef},$cutPos;

		while($pos>$outerBinEnd){ # filled outer region, start kernel density analysis
			$self->cleanUp($arrInnerRef,$innerBinStart); # remove reads of last step (window)
			$self->cleanUp($arrOuterRef,$outerBinStart); # remove reads of last step (window)

			my $outerReadNum=@{$arrOuterRef}-1; # get number of reads of outer bin
			my $innerReadNum=0; # number of reads of inner bin
#			my $temNum=@{$arrInnerRef}-1;

			for(my $x=0;$x<@{$arrInnerRef};$x++){ # this code is to solve the problem that reads is sorted by the left most position rather thant the 5' end position
				last if $$arrInnerRef[$x]-100>$innerBinEnd;
				$innerReadNum++ if $$arrInnerRef[$x]>$innerBinStart && $$arrInnerRef[$x]<$innerBinEnd;
			}

			my $mappable=0;
			my $innerMappable;
			my $normReadNum=0;

			if($innerReadNum>0){
				($mappable,$innerMappable)=$self->overlap(\@mappedSt,\@mappedEd,$outerBinStart,$outerBinEnd,$innerBinStart,$innerBinEnd); # SLOW
				if($mappable>0){
					my $val=$innerReadNum-int($innerBinSize*$outerReadNum/$mappable);
					if($val>0){
						$normReadNum=$val;
					}
				}
			}

#			print BED "$chr\t$pos\t",$pos+1,"\t$innerReadNum\_$normReadNum\_$mappable\_$outerReadNum\_$innerBinSize\n";
			$normReadNum=0 if $normReadNum<0;


			if($normReadNum>0){
				for(my $x=0;$x<$innerBinSize;$x++){ # gaussion smooth
				$$densityRef[$x]+=$normReadNum*$$FixedArr[$x];;
	     	   }
			}

			my @stepDensity=splice(@{$densityRef}, 0, $stepSize); # density in the step (window)

			# get peak
			if($normReadNum>=$cutoff){
				$peakTh=1;
			}

			my $sum=0;
			for(my $x=0;$x<@stepDensity;$x++){
#				print "$chr\t",$innerBinStart+$x,"\t",$innerBinStart+$x+1,"\t$stepDensity[$x]\_$normReadNum\n";
				$sum+=$stepDensity[$x];
				if($stepDensity[$x]>=$bound){ # threshold should be determined before and sd
					if($peakOn==0){ # peak initiate
						$peakOn=1;
						$peakStart=$innerBinStart+$x;
						$peakEnd=$innerBinStart+$x;
						$peakScore=$normReadNum;
					}
					else{ # peak extend
						$peakEnd++;
						$peakScore=$normReadNum if $peakScore<$normReadNum;
					}
				}
				else{
					if($peakOn==1){ # peak end
						my $pv=$$pval{$peakScore} || 0;
						print BED "$chr\t$peakStart\t$peakEnd\t$peakScore\tpeak\t$pv\n" if $peakTh==1;
						$peakOn=0;
						$peakTh=0;
					}
				}
				$stepDensity[$x]+=0.00000001; # for comparible to pybigwig
			}

			# get wig out
#            if($$arrInnerRef[0]>$innerBinEnd){# no reads, make a gap so that the data of this gap will not be dumped in the wig file

			if($sum==0){
				if($gap==0){ # gap start and end of current entry, output values if any, then defined the start point of next entry
#					print "$dbg we got 2 $bwStart, @bwEntries\n\n";
					if(@bwEntries>0){
#						print "adding $Nchr, $bwStart, @bwEntries\n";
#						print "we got 2 $bwStart, \n\n";
						$bw->addEnt($Nchr,$bwStart,\@bwEntries,1,1);
					}
					$gap=1;
					undef @bwEntries;
					$bwStart=$innerBinStart;
				}
				else{ # gap extend, defined the start point of next entry
					$bwStart=$innerBinStart;
				}
			}
			else{
				if($gap==0){ # extend current entry
					push @bwEntries,@stepDensity;
				}
				else{ # gap end, start new entry
					push @bwEntries,@stepDensity;
					$bwStart=$innerBinStart;
					$gap=0;
				}
			}


			$innerBinStart+=$stepSize;
			$innerBinEnd+=$stepSize;
			$outerBinStart+=$stepSize;
			$outerBinEnd+=$stepSize;
		}
#		print "output the loop\n";
	}

    if(@bwEntries>0){
#		print "we got 1 $bwStart, @bwEntries\n";
        $bw->addEnt($Nchr,$bwStart,\@bwEntries,1,1);
        undef @bwEntries;
    }

	$bw->close() if $aa>1;
	close BED;

}

sub cleanUp{
    my $self=shift;
    my $arr=shift;
    my $st=shift;
    while($$arr[0]<$st){
		last if !$$arr[0];
        shift @{$arr};
    }
}


sub overlap{
    my $self=shift;

    my $arrSt=shift;
    my $arrEd=shift;
    my $st=shift;
    my $ed=shift;
	my $inSt=shift;
	my $inEd=shift;
    my $sum=0;
	my $inSum=0;
    my $x=0;

    while(1){
        last if !$$arrEd[$x];
        if($$arrEd[$x]<$st){
            shift @{$arrSt};
            shift @{$arrEd};
        }
        elsif($$arrSt[$x]>$ed){
            last;
        }
        else{
            my $up=$$arrSt[$x]>$st ? $$arrSt[$x] : $st;
            my $down=$$arrEd[$x]<$ed ? $$arrEd[$x] : $ed;
            $sum+=$down-$up+1;
            $x++;
        }
    }
    return ($sum, $inSum);
}


sub getmapReg{
    my $self=shift;

    my $file=shift;
    my $hash;
    open F,"less $file | " or die ("cannot open $file\n");
    while(<F>){
        my @tem=split;
        $$hash{$tem[0]}{$tem[1]}=$tem[2];
    }
    close F;

    return $hash;
}

sub timer {
	my $pos=shift;
	my $aa=shift;
	my $tt=shift;
	my $genomeSize=shift;
	$aa++;
	print STDERR "got $aa reads done\n" if !($aa%1000000);
	if($pos>1000000*$aa){
		$aa++;
		my $ntt=time;
		my $utt=$ntt-$tt;
		$tt=$ntt;

#		my $rest=($genomeSize-$pos)*$utt/3600/1000000;
#		print STDERR " $utt second per Mb, need $rest hours for $genomeSize bp\n";
#		print STDERR "$ntt\n";
	}
	return ($tt,$aa);
}
use Inline Python => <<'END';


import pyBigWig

class bwObj:

    def __init__(self):
#		self.start="new bw object being created"
        print "new bw object being created"

    def addHead(self,chr="10",chrLen=150982314):
        print "PY add header", chr, chrLen
        self.bw.addHeader([(str(chr),int(chrLen))])

    def addEnt(self,chr="10",st=900,values=[0,0,0],sp=20,ste=30):
#       print "add entries", chr, st, values, sp, ste
        self.bw.addEntries(str(chr), st, values=values, span=sp, step=ste)

    def openBw(self,file,mode="r"):
        self.bw=pyBigWig.open(file,mode)
#       print "open file"
        return self.bw;

    def close(self):
        self.bw.close()

END



1;
