#! /usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

my ($file);
my ($outdir,$prefix) = ("./","All");

GetOptions (
	"f:s"	=>	\$file,
	"p:s"	=>	\$prefix,
	"o:s"	=>	\$outdir
);

sub usage {
	die "
 perl $0

 -f [s]	All Samples MetaPhlAn2 Merge file
 -p [s]	Prefix of output [default All]
 -o [s]	Output directory [default .]\n\n";
}

unless (defined $file)
{
	&usage;
	exit;
}

chomp (my $pwd = `pwd`);
$outdir = ($outdir eq "." || $outdir eq "./") ? $pwd : ($outdir =~ /^\//) ? $outdir : $pwd."/".$outdir;
system ("mkdir -m 755 -p $outdir") unless (-d $outdir);

my (@sample,%stat,%total);
my @level = ("Kingdom","Phylum","Class","Order","Family","Genus","Species");
chomp(my $sum=`grep -v "|" $file |awk -F '\t' '{ sum += \$2;}; END { print sum }'`);
open ABU,">$file.reAbu.xls";
open FILE, "$file" || die "can not open: $file, $!\n";
my $id=(split /\//,$file)[-1];
$id=~s/.Kraken2.out//;
print ABU "#SampleID\t$id\n";
push @sample,"head",$id;
while (<FILE>)
{
	chomp;
	if (/^#/)
	{
	} else {
		next unless (/^d__/);
		my @array = split /\t+/;
		my $reAbu=$array[1]/$sum*100;
		$array[0]=~s/ /_/g;
		$array[0]=~s/_sp\._/_sp_/g;
		print ABU "$array[0]\t$reAbu\t$array[1]\n";
		my @species = split /\|/,$array[0];
		for (my $i=1; $i<=$#array; $i++)
		{
			$stat{$level[0]}{$sample[$i]} ++ if (@species == 1 && $array[$i] > 0);
			$stat{$level[1]}{$sample[$i]} ++ if (@species == 2 && $array[$i] > 0 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
			$stat{$level[2]}{$sample[$i]} ++ if (@species == 3 && $array[$i] > 0 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
			$stat{$level[3]}{$sample[$i]} ++ if (@species == 4 && $array[$i] > 0 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
			$stat{$level[4]}{$sample[$i]} ++ if (@species == 5 && $array[$i] > 0 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
			$stat{$level[5]}{$sample[$i]} ++ if (@species == 6 && $array[$i] > 0 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
			$stat{$level[6]}{$sample[$i]} ++ if (@species == 7 && $array[$i] > 0 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
		}
		$total{$level[0]} ++ if (@species == 1);
		$total{$level[1]} ++ if (@species == 2 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
		$total{$level[2]} ++ if (@species == 3 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
		$total{$level[3]} ++ if (@species == 4 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
		$total{$level[4]} ++ if (@species == 5 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
		$total{$level[5]} ++ if (@species == 6 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
		$total{$level[6]} ++ if (@species == 7 && $species[-1] !~ /noname$/ && $species[-1] !~ /unclassified$/);
	}
}
close FILE;
close ABU;
#print Dumper (\%stat);

open OUT, ">$outdir/$prefix\_Species_stat.xls";
print OUT "#Single Sample Statistics Result:\n";
shift (@sample);
my $title;
foreach my $sample (@sample)
{
	$title .= "\t$sample";
}
print OUT "#SampleID"."$title\n";

foreach my $level (@level)
{
	my $print = "";
	foreach my $sample (@sample)
	{
		
		my $level_num = (defined $stat{$level}{$sample}) ? $stat{$level}{$sample} : 0;
		$print .= "\t$level_num";
	}
	print OUT "$level"."$print\n";
}

if (@sample >= 2)
{
	print OUT "\n#All Samples Statistics Result:\n";
	foreach my $level (@level)
	{
		print OUT "$level\t$total{$level}\n";
	}
}
close OUT;
