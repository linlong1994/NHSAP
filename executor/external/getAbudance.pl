#!/usr/bin/perl -w
use strict;
use File::Basename;
die "perl $0 bamstat read_classification.txt reference_fai reference_name\n" unless @ARGV==4;


my %map;
my $stat=shift;
open Stat,"less $stat|";
while(<Stat>){
	chomp;
	my @t=split /\t/;
	$map{$t[0]}=$t[3];
}
close Stat;

open Read,shift;
my $Refai = shift;
my $reference_len = `awk '{sum+=\$2} END {print sum}' $Refai`;
chomp $reference_len;
my $refName=shift;
while(<Read>){
	chomp;
	my @t=split /\t/;
	my $count=0;
	for my $j(1..$#t){
	$count +=$1 if $t[$j]=~/\:(\d+)$/;
	}
	my $abu;
	if ((! exists($map{$t[0]})) or $map{$t[0]} eq 0 ){
		$abu = 'NA'
	}
	else{
	$abu=2*$count*3099706404/($reference_len*$map{$t[0]});}
	print "$t[0]\t$refName\t$abu\n";
}
close Read;
