#!/usr/bin/perl -w
use strict;
use Getopt::Std;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

$,= " ";

my %OPTS;
getopts('i:j:k:l:m:n:o:s:a:',\%OPTS);

my $out = $OPTS{o};

my $fileOut1=$out.".mutation_assignment.txt";
my $fileOut2=$out.".subclonal_structure.txt";

my $summaryI=$OPTS{i};
my $snvFileI=$OPTS{j};

my $summaryJ=$OPTS{k};
my $snvFileJ=$OPTS{l};

my $summaryK=$OPTS{m};
my $snvFileK=$OPTS{n};

my $headerRef;
my $sumRef;

($headerRef,$sumRef)=readSummary($summaryI);
my @headerI=@$headerRef;
my %sumI=%$sumRef;

($headerRef,$sumRef)=readSummary($summaryJ);
my @headerJ=@$headerRef;
my %sumJ=%$sumRef;

($headerRef,$sumRef)=readSummary($summaryK);
my @headerK=@$headerRef;
my %sumK=%$sumRef;

printSummary(\@headerI,\%sumI);
printSummary(\@headerJ,\%sumJ);
printSummary(\@headerK,\%sumK);

## get delta log likelihood;

my $line="# n cna-llh baf-llh snv-llh total-llh total-bic";

my @scoreI=split(/\s+/,$sumI{$line}[0]);
my @scoreJ=split(/\s+/,$sumJ{$line}[0]);
my @scoreK=split(/\s+/,$sumK{$line}[0]);

my $dL1=$scoreJ[-3]-$scoreI[-3];
my $dL2=$scoreK[-3]-$scoreJ[-3];

print $dL1, $dL2,"\n";

# cut = total amount of high confidence posterior prob (0.9) needed to justify a cluster 
# cutL = Log-likelihood based model selection
 
my $cut=$OPTS{a}; 
my $cutL=$OPTS{s};

print "Parameters", $cut,$cutL,"\n";

my $model;
my $cloneSizes;
if($dL1> $cutL ){
	if($dL2> $cutL){
		($headerRef,$sumRef)=readSNVs($snvFileK);
		$model="# 3 clones";
		$cloneSizes=$sumK{$model}[0];
	}
	else{
		($headerRef,$sumRef)=readSNVs($snvFileJ);
		$model="# 2 clones";
		$cloneSizes=$sumJ{$model}[0];
	}
}
else{
	($headerRef,$sumRef)=readSNVs($snvFileI);
	$model="# 1 clones";
	$cloneSizes=$sumI{$model}[0];
}

my @header=@$headerRef;
my %values=%$sumRef;

hqStates(\%values,$cut,$cloneSizes);

sub hqStates{
	# count how much posterior is in a well supported state for mutations where a clone is in zero state
	my ($sumRef,$cut,$cloneSizes)=@_;
	my %values=%$sumRef;
	my @states=split(/\s+/,$values{"#copynumbers:"}[0]);
	my @statesCG;
	my @freqs=split(/\s+/,$values{"#clonal frequencies:"}[0]);
	my @snv=@{$values{"#chr locus PostDist"}};
	my @sumP=();
	
	my %clusters;
	my %clustersCG;
	
	for(my $i=0; $i<scalar(@states); $i++){
		push(@sumP,0.0);
	}
	for(my $i=0; $i<scalar(@states); $i++){
		my @genotype=split(//,$states[$i]);
		my $gt="";
		for(my $j=0; $j<scalar(@genotype); $j++){
			if($genotype[$j]>0){
				$gt=$gt."1";
			}
			else{
				$gt=$gt."0";
			}
		}
		push(@statesCG,$gt);
	}
	
	# iterate over SNVs and integrate    
	for(my $l=0; $l<scalar(@snv); $l++){
		my @fields=split(/\s+/,$snv[$l]);
		
		for(my $i=0; $i<scalar(@states); $i++){
			$sumP[$i]+=$fields[2+$i];
			if(defined $clusters{$statesCG[$i]}){
				$clusters{$statesCG[$i]}+=$fields[2+$i];
				if($fields[2+$i]>0.9){
					$clustersCG{$statesCG[$i]}+=$fields[2+$i];
				}
			}
			else{
				$clusters{$statesCG[$i]}=$fields[2+$i];
				if($fields[2+$i]>0.9){
					$clustersCG{$statesCG[$i]}=$fields[2+$i];
				}
			}
		}
	
	}
	
	my @cloneM=();
	my $nClones=scalar(split(//,$states[0]));
	my @cloneSizesV=split(/\s/,$cloneSizes);
	my %clusterFreq;
	
	while( my( $key, $value ) = each %clusters ){
		my @gt=split(//,$key);
		my $freq=0.0;
		for(my $i=0; $i<scalar(@gt); $i++){
			if($gt[$i]==1){
				$freq+=$cloneSizesV[$i];
			}
		}
		$clusterFreq{$key}=$freq;
	};
	
	my @bestI=();
	foreach (sort {$clusterFreq{$b} <=> $clusterFreq{$a}} keys %clusterFreq) {
		# do not print out zero cluster
		if($clustersCG{$_}>$cut){
			print "$cut Included $_: $clusterFreq{$_} $clusters{$_} $clustersCG{$_}\n";
			push(@bestI,$_);
		}
	}
	
	while( my( $key, $value ) = each %clusters){
		if(defined $clustersCG{$key}){
			print "$key $value $clustersCG{$key} $clusterFreq{$key}\n";
		}
		else{
			print "$key $value 0.0 $clusterFreq{$key}\n";
		}
	}
	
	$,= "\t";
	open(my $fh, '>', $fileOut1) or die "Could not open file '$fileOut1' $!";
	print $fh "#\t\t";
	for(my $i=0; $i<scalar(@bestI); $i++){
		print $fh $bestI[$i]." ".$clusterFreq{$bestI[$i]},"";
	}
	print $fh "\n";
	
	print $fh "#chr\tlocus\t";
	for(my $i=0; $i<scalar(@bestI); $i++){
		print $fh "cluster".($i+1),"";
	}
	print $fh "\n";
	
	for(my $l=0; $l<scalar(@snv); $l++){
		my @fields=split(/\s+/,$snv[$l]);
		
		my %snvCG;
		
		# rescale probabilities to take into accout removal of CG clusters
		my $sum=0.0;
		for(my $i=0; $i<scalar(@states); $i++){
			
			if(defined $snvCG{$statesCG[$i]}){
				$snvCG{$statesCG[$i]}+=$fields[2+$i];
				$sum+=$fields[2+$i];
			}
			else{
				$snvCG{$statesCG[$i]}=$fields[2+$i];
				$sum+=$fields[2+$i];
			}
		}
		
		print $fh $fields[0],$fields[1],"";
		for(my $i=0; $i<scalar(@bestI); $i++){
			if($sum!=0.0){
				print $fh $snvCG{$bestI[$i]}/$sum,"";
			}
			else{
				print $fh 1.0/scalar(@bestI),"";
			}
		}
		print $fh "\n";
	}
	close($fh);
	
	open($fh, '>', $fileOut2) or die "Could not open file '$fileOut2' $!";
	print $fh "#cluster \t n_ssms \t proportion \t\n";
	for(my $i=0; $i<scalar(@bestI); $i++){
		print $fh ($i+1),$clusters{$bestI[$i]},$clusterFreq{$bestI[$i]};
		print $fh "\n";
	}
	close($fh);
    
}

sub readSNVs{
	my ($file) = @_;
	my @header;
	my %values;
	
	open my $IN_FILE, '<', $file or die "Cannot open  input file";
	foreach my $line (<$IN_FILE>) {
		chomp($line);
		
		my @fields=split(/\s+/,$line);
		if(($line =~ /\#/)){
			push(@header,$line);
		}
		else{
			push(@{$values{$header[-1]}},$line);
		}
	};
	close $IN_FILE;
	return (\@header,\%values);
}

sub printSummary{
	my ($ref1,$ref2) = @_;
	my @headerI=@$ref1;
	my %sumI=%$ref2;
	
	for(my $i=0; $i<scalar(@headerI); $i++){
		for(my $j=0; $j<scalar(@{$sumI{$headerI[$i]}}); $j++){
			print $headerI[$i],$sumI{$headerI[$i]}[$j],"\n"; 
		}
	}
}

sub readSummary {
	my ($file) = @_;
	
	my @headerI;
	my @headerJ;
	
	my %sumI;
	my %sumJ;
	
	open my $IN_FILE, '<', $file or die "Cannot open input file";
	foreach my $line (<$IN_FILE>) {
		chomp($line);
		
		my @fields=split(/\s+/,$line);
		if(($line =~ /\#/)){
			push(@headerI,$line);
		}
		else{
			push(@{$sumI{$headerI[-1]}},$line);
		}
	};
	close $IN_FILE;
	return (\@headerI,\%sumI); 
};