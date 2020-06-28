
sub unmapped_fq_from_bam{

    my $samtools = $_[0];
    my $bam = $_[1];
    my $sample_id = $_[2];
    my $outdir = $_[3];

    open(I,"$samtools view -f 4 $bam|") or die "can not open bamfile"; 
    open OUT,"| gzip > $outdir/$sample_id.fq.gz" or die "can not write to gz";
    while(<I>){
        chomp;
        @samcolumns=split;
        if($samcolumns[2]  eq "*" ){
            print OUT "@".$samcolumns[0]."\n";
            print OUT "$samcolumns[9]\n";
            print OUT "+\n";
            print OUT "$samcolumns[10]\n";
        }
    } 
    close(I);
    close(OUT);
}


$samtools = $ARGV[0];
$bam = $ARGV[1];
$sample_id = $ARGV[2];
$outdir = $ARGV[3];


unmapped_fq_from_bam($samtools,$bam,$sample_id,$outdir);