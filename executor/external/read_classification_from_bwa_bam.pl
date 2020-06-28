
sub read_classification_from_bwa_bam{

    my $samtools = $_[0];
    my $bam = $_[1];
    my $sample_id = $_[2];
    my $outdir = $_[3];

    open OU,">$outdir/$sample_id.read_classification.txt" or die $!;
    open(I,"$samtools view $bam | awk '\$3!=\"*\" {print \$0}' |");
    my %reads_ref=();
    my %reads_evalue=();
    while($fileline=<I>){
        next if(/^\s$/);
        @line=split(/\s+/,$fileline);
        unless(exists $reads_ref{$line[0]}){
            $reads_ref{$line[0]}=$line[2];
            $reads_evalue{$line[0]}=$line[4]
        }
        else{
            if($line[4]>$reads_evalue{$line[0]}){
                $reads_ref{$line[0]}=$line[2];$reads_evalue{$line[0]}=$line[4]
            }
        }
    }
    my %read_classification_persample=();
    while(($reads,$ref)=each %reads_ref){
        unless(exists $read_classification_persample{$ref}){
            $read_classification_persample{$ref}=1
        }
        else{
            $read_classification_persample{$ref}+=1
        }
    }
    print OU "$sample_id";
    foreach $ref (sort { $read_classification_persample{$b} <=> $read_classification_persample{$a} } keys  %read_classification_persample){
        print OU "\t$ref:$read_classification_persample{$ref}";
    }
    print OU "\n";
    close OU;
}


$samtools = $ARGV[0];
$bam = $ARGV[1];
$sample_id = $ARGV[2];
$outdir = $ARGV[3];

read_classification_from_bwa_bam($samtools,$bam,$sample_id,$outdir);
