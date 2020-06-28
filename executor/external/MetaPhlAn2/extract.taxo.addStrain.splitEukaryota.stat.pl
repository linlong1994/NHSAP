#!/usr/bin/perl -w
use strict;
die "perl $0 /hwfssz5/ST_MCHRI/BIGDATA/USER/luoguangwen/bin/NIPT/metaphlan2/fungi.cut2.xls /hwfssz5/ST_MCHRI/BIGDATA/USER/luoguangwen/bin/NIPT/metaphlan2/parasite.cut2.xls MetaPhlAn2.out outdir > MetaPhlAn2.out.order\n" unless @ARGV;
my %fungi;
my %para;
my %h;
my %h0;
my %spe;
my %st;
open Fungi,shift;
while(<Fungi>){
	chomp;
	$fungi{$_}=1;
}
close Fungi;
open Para,shift;
while(<Para>){
	chomp;
	$para{$_}=1;
}
close Para;

my %wh;
open F,shift;
chomp(my $dir=`pwd`);
$dir=shift if @ARGV;
my $H=<F>;
$H=~s/\.MetaPhlAn2//g;
$H=~s/\#SampleID//;
$H=~s/^ID//;
print $H;
my $euk=0;
my $Para=0;
my $Fungi=0;
my $opt=0;
my ($FungiCOV,$FungiAVEL,$FungiReads)=(0,0,0);
my ($ParasiteCOV,$ParasiteAVEL,$ParasiteReads)=(0,0,0);
while(<F>){
	next if /^\#/;
	my @t=split /\t/,$_,2;
	my @s=split /\|/,$t[0];
	my $Kingdom=$s[0]; #add
	$Kingdom=~s/^k__//; #add
	my $c=scalar(@s);
	my $id=$s[-1];
	$id=~s/^k__|^p__|^c__|^o__|^f__|^g__|^s__|^t__//g;
	$h0{$c} .=$_;
	$h{$c} .="$id\t$t[1]" unless($c==1 && $id eq "Eukaryota");
	$euk =1 if($c==1 && $id eq "Eukaryota");
	if($euk ==1 && $c==7){
		my @ss=split /\t/,$t[1];
		my($abu,$coverage,$average_genome_length_in_the_clade,$estimated_number_of_reads_from_the_clade);
		if(@ss==4){
			$opt=1;
			($abu,$coverage,$average_genome_length_in_the_clade,$estimated_number_of_reads_from_the_clade)=@ss[0,1,2,3];
		}else{
			$abu=$ss[0];
		}
		if (exists $fungi{$id}){
			$Fungi+=$abu;
			$Kingdom="Fungi";
			if($opt==1){
				$FungiCOV+=$coverage;
				$FungiAVEL+=$average_genome_length_in_the_clade;
				$FungiReads+=$estimated_number_of_reads_from_the_clade;
			}
		}elsif (exists $para{$id}){
			$Para +=$abu;
			$Kingdom="Parasite";
			if($opt==1){
				$ParasiteCOV+=$coverage;
				$ParasiteAVEL+=$average_genome_length_in_the_clade;
				$ParasiteReads+=$estimated_number_of_reads_from_the_clade;
			}
		}
		$spe{$Kingdom}++; #add
	}elsif($c==7){
		$spe{$Kingdom}++;#add
	}
	$st{$c} ++;
	$wh{$c} .="$Kingdom\|$id\t$t[1]" if $c==7; #add
	if($c==3 && $id=~/unclassified$/){
		$h{4} .="c__$id\t$t[1]";
		$h{5} .="c__$id\t$t[1]";
		$h{6} .="c__$id\t$t[1]";
		$h{7} .="c__$id\t$t[1]";
		$h{8} .="c__$id\t$t[1]";
		$st{4} ++;
		$st{5} ++;
		$st{6} ++;
		$st{7} ++;
		$st{8} ++;
		$wh{7} .="$Kingdom\|c__$id\t$t[1]";#add
		$spe{$Kingdom}++;#add
	}
	if($c==4 && $id=~/unclassified$/){
		$h{5} .="o__$id\t$t[1]";
		$h{6} .="o__$id\t$t[1]";
		$h{7} .="o__$id\t$t[1]";
		$h{8} .="o__$id\t$t[1]";
		$st{5} ++;
		$st{6} ++;
		$st{7} ++;
		$st{8} ++;
		$wh{7} .="$Kingdom\|o__$id\t$t[1]";#add
		$spe{$Kingdom}++;#add
	}
	if($c==5 && $id=~/unclassified$/){
		$h{6} .="f__$id\t$t[1]";
		$h{7} .="f__$id\t$t[1]";
		$h{8} .="f__$id\t$t[1]";
		$st{6} ++;
		$st{7} ++;
		$st{8} ++;
		$wh{7} .="$Kingdom\|f__$id\t$t[1]";#add
		$spe{$Kingdom}++;#add
	}
	if($c==6 && $id=~/unclassified$/){
		$h{7} .="g__$id\t$t[1]";
		$h{8} .="g__$id\t$t[1]";
		$st{7} ++;
		$st{8} ++;
		$wh{7} .="$Kingdom\|g__$id\t$t[1]";#add
		$spe{$Kingdom}++;#add
	}
	if($c==7 && $id=~/unclassified$/){
		$h{8} .="s__$id\t$t[1]";
		$st{8} ++;
	}
}
open K,"> $dir/Kingdom.rabun.xls";
open P,"> $dir/Phylum.rabun.xls";
open C,"> $dir/Class.rabun.xls";
open O,"> $dir/Order.rabun.xls";
open F,"> $dir/Family.rabun.xls";
open G,"> $dir/Genus.rabun.xls";
open S,"> $dir/Species.rabun.xls";
open T,"> $dir/Strain.rabun.xls";
print K $H;
print P $H;
print C $H;
print O $H;
print F $H;
print G $H;
print S $H;
print T $H;

my $tmp1="";
$tmp1="\t$FungiCOV\t$FungiAVEL\t$FungiReads" if $opt==1;
my $tmp2="";
$tmp2="\t$ParasiteCOV\t$ParasiteAVEL\t$ParasiteReads" if $opt==1;
$h{1} .="Fungi\t$Fungi"."$tmp1\n" if $Fungi>0;
$h{1} .="Parasite\t$Para"."$tmp2\n" if $Para>0;
for (sort keys %h){
	print $h0{$_};
	print K $h{$_} if $_==1;
	print P $h{$_} if $_==2;
	print C $h{$_} if $_==3;
	print O $h{$_} if $_==4;
	print F $h{$_} if $_==5;
	print G $h{$_} if $_==6;
	print S $h{$_} if $_==7;
	print T $h{$_} if $_==8;
}
open ST,"> $dir/Taxo.stat.txt";
print ST "Kingdom\t$st{1}
Phylum\t$st{2}
Class\t$st{3}
Order\t$st{4}
Family\t$st{5}
Genus\t$st{6}
Species\t$st{7}
Strain\t$st{8}\n";

open WS,"> $dir/Species.rabun.addKingdom.xls";
print WS $H;
print WS $wh{7};
open STT,"> $dir/species.stat.txt";
for (sort keys %spe){
	print STT "$_\t$spe{$_}\n";
}
