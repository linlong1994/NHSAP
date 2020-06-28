#!/usr/bin/perl
die "perl $0 Kraken2.out.reAbu.xls > Kraken2.Species.reAbu.xls 2> Kraken2.Species.addKingdom.reAbu.xls\n" unless @ARGV;
my %h;
my %h0;
open F,shift;
my $H=<F>;
$H=~s/\#SampleID//;
print $H;
print STDERR $H;
while(<F>){
	my @t=split /\t/,$_,2;
	my @s=split /\|/,$t[0];
	my $id=$s[-1];
	next unless $id=~/^s__/;
	$id=~s/^s__//g;
	$s[0]=~s/^d__//g;
	print "$id\t$t[1]";
	print STDERR "$s[0].$id\t$t[1]";
}
close F;
