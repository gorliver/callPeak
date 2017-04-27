package threshold;
use strict;
use warnings;

#use parent qw(smooth);

use Bio::DB::HTS::Tabix;

sub new{
	my $class=shift;
	my $self={};
	return bless $self,$class;
}

sub getData{
	my $self=shift;
	my $name=shift;
	if($self->{$name}){
		return $self->{$name};
	}
	else{
		die("no data for $name\n");
	}
}

sub setChromInfo{
	my $self=shift;
	my $chromFile=shift;
	open my $cf,$chromFile or die ("no chromosome size file\n");
	my $hash;
	my $idx=0;
	my $genomeSize=0;
	while(<$cf>){
		$idx++;
		my @tem=split;
		$$hash{$idx}=\@tem;
		$genomeSize+=$tem[1];
	}
	close $cf;
	$self->{'ChromInfo'}=$hash;
	$self->{'GenomeSize'}=$genomeSize;
}


# random get 500 regions from 1 to 
sub makeRandPos{
	my $self=shift;
	my $randNum=shift;
	my $halfIn=shift;
	my $halfOut=shift;
	my $genomeSize=$self->getData('GenomeSize');
	my $chrLenHash=$self->getData('ChromInfo');

	my $innerRead;
	my $outerRead;
	
	my %randPos;
	srand 1;
	my $cnt=0;
	print STDERR "we got genome size if $genomeSize\n";
	for(my $x=0;$x<$randNum;$x++){
		my $rp=int(rand($genomeSize));
		next if $randPos{$rp};
		$randPos{$rp}++;

		my ($asChr,$asPos)=$self->assign($rp);
		next if $asPos<$halfOut;
		$cnt++;
	
		$$innerRead{$asChr}{$asPos}=0;
		$$outerRead{$asChr}{$asPos}=0;
	}
	
	$self->{'InnerRead'}=$innerRead;
	$self->{'OuterRead'}=$outerRead;
}

sub assign{
	my $self=shift;
    my $num=shift;
	my $chrLenRef=$self->getData('ChromInfo');

	my $stChrIdx=1;
	my $stChrName=$$chrLenRef{$stChrIdx}->[0];
	my $stChrLen=$$chrLenRef{$stChrIdx}->[1];
	my $stPos=1;
	my $edPos=$stPos+$stChrLen-1;
	while(1){
		if($num<$edPos){
			my $assignChr=$stChrName;
			my $assignPos=$num-$stPos;
			return ($assignChr,$assignPos);
		}
		if($num>$edPos){
			$stPos=$edPos;
			$stChrIdx++;
			$edPos+=$$chrLenRef{$stChrIdx}->[1];
			$stChrName=$$chrLenRef{$stChrIdx}->[0];
			$stChrLen=$$chrLenRef{$stChrIdx}->[1];
		}
	}
}
	
## determine the distribution of number of reads in each window

sub makeDistr{	
	my $self=shift;
	my $bamFile=shift;
	
	my $halfIn=shift;
	my $halfOut=shift; #shift

	my $uniqRegionFile=shift;
		
	my $innerRead=$self->getData('InnerRead');
	my $outerRead=$self->getData('OuterRead');
		
	#print "#rand is @sortRandPos\n";
	
	my $Nchr=1000;
	
	my $stPoint=0;
	my $next=$stPoint;
	
	my @sortRandPos; # =@{$self->getData('RandPos')};

	my $aa=0;	
#	open my $bamFile,$ARGV[0] or die;


	print STDERR "read bam\n";
	open my $bamHandle,"samtools view -h $bamFile | " or die ("cannot open $bamFile for makeDistr\n");
	while(<$bamHandle>){
		next if /^\@/;
		$aa++;
#		my ($chr,$pos)=(split)[2,3];
		my ($chr,$flag,$pos,$qual,$tag)=(split)[2,1,3,4,5];
		next if $qual<20;
		next if $chr=~/sca/;
		next if $chr eq "*";

	    if($chr ne $Nchr){#reach new chr or start from the first chr, initiate the rand pos array;
	        $Nchr=$chr;
			@sortRandPos=sort {$a<=>$b} keys %{$innerRead->{$Nchr}};
			my $nn=@sortRandPos;
			$stPoint=0;
			$next=0;
#			print "@sortRandPos\n";
			print STDERR "#### start $chr for $nn pos\n";
	        redo;
		}
		
		next if !$sortRandPos[$stPoint];
		my $rdPos=$sortRandPos[$stPoint];
		my $innerSt=$rdPos-$halfIn;
		my $innerEd=$rdPos+$halfIn;
		my $outerSt=$rdPos-$halfOut;
		my $outerEd=$rdPos+$halfOut;
	
		if(!($aa%1000000)){
			print STDERR "finished $aa reads\n";
#			print STDERR time(),"\n";
		}
	
		if($pos>$outerEd){
			$stPoint++;
			$next=$stPoint;

			redo;
		}
		my $next=$stPoint;
#		if($stPoint>12670){
#			print  "$aa,$chr,$Nchr,$aa $pos \t $rdPos\t $outerSt $outerEd\tstart is $stPoint, $next\n";
#		}


		while($pos<$outerEd){
#			print "\tstart iterative($next) $aa $Nchr $pos $outerSt $outerEd\n" if $Nchr==2;
			
			if($pos<$outerSt){
				last;
			}
#			print "\tadd one\n";
			$$outerRead{$Nchr}{$rdPos}++;
			$$innerRead{$Nchr}{$rdPos}++ if $pos>=$innerSt && $pos<$innerEd;

			$next++;
#			print "\t next if $next, iterative\n";
	        last if !$sortRandPos[$next];
	        $rdPos=$sortRandPos[$next];
	        $innerSt=$rdPos-$halfIn;
	        $innerEd=$rdPos+$halfIn;
	        $outerSt=$rdPos-$halfOut;
	        $outerEd=$rdPos+$halfOut;
		}
	}
	close $bamHandle;
	my $tabix=Bio::DB::HTS::Tabix->new(filename=>$uniqRegionFile);
#	print "#output\n";
	my @density;
#	foreach my $one(@sortRandPos){
	foreach my $oneChr(keys %{$innerRead}){
		foreach my $onePos(keys %{$innerRead->{$oneChr}}){
#		print"we got $one, $innerRead{$one}, $outerRead{$one}, $halfIn, $halfOut\n";

			## get mapped regions
			my $mappable=mapReg($tabix,$oneChr,$onePos-$halfOut,$onePos+$halfOut);
			my $normRead=$$innerRead{$oneChr}{$onePos}-int(2*$halfIn*$$outerRead{$oneChr}{$onePos}/2/$halfOut);
#			my $normRead=$innerRead{$one}-int(2*$halfIn*$$outerRead{$one}/2/$halfOut);
			$normRead=0 if $normRead<0;
#		print "#output is $oneChr\t$onePos\t$$innerRead{$oneChr}{$onePos}\t$$outerRead{$oneChr}{$onePos}\t$normRead\t$mappable\n";
		
			push @density,$normRead;
		}
	}
	@density=sort {$a<=>$b} @density;
	$self->{'Density'}=\@density;
}

#sub getDistr{
#	my $self=shift;
#	return $self->{'Density'}
#}
## kernel density smooth and determine the threshold

sub mapReg{
	my $tab=shift;
	my $chr=shift;
	my $st=shift;
	my $ed=shift;
	my $reg="$chr:$st-$ed";
	my $iter=$tab->query($reg);
	my $nn;
	while(my $l=$iter->next){
		my @tt=split /\t/,$l;
		my $up=$tt[1]>$st ? $tt[1] : $st ;
		my $down=$tt[2]<$ed ? $tt[2] : $ed;
		$nn+=$down-$up;
	}
	return $nn;
}

sub calPval{
	my $self=shift;
	my $sumx;
	my $sumx2;
	my @density=@{$self->getData('Density')};
	my $pvalue;
	
	
	my $count=@density;
	foreach my $one(@density){
		$sumx+=$one;
		$sumx2+=$one*$one;
	}
	my $sig=sqrt($sumx2/$count-($sumx/$count)**2);
	my $w=1.06*$sig*($count**(-1/5));
#	print "w is $w\n";
	
	for(my $x=0;$x<$density[-1];$x+=1){
		my $prob=0;
		foreach my $val(@density){
			my $z=$x-$val;
			my $pp=-0.5*(($z*$z/$w)**2);
			my $ee=exp($pp);
			my  $ff=$ee/$w/sqrt(2.0*3.1415926);
			$prob+=$ff;
		}
		$prob=$prob/$count;
#		print "$x\t$prob\n";
		$$pvalue{$x}=$prob;
		last if $prob<1e-20;
	}
	$self->{'Pvalue'}=$pvalue;
}

1;

