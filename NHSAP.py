
import time 
import sys
import argparse
import ast

def gnhs(args):
    from executor.executor import GnhsRunner

    gsr = GnhsRunner(args)

    result = gsr.gnhs_runner()

    is_success = False
    if result:
        is_success = True
    # for p in processer:
    #     if p.exitcode != 0:
    #         is_success = False

    return is_success

def gnha(args):
    from executor.executor import GnhaRunner

    gar = GnhaRunner(args)

    processer = gar.gnha_runner()

    is_success = True
    for p in processer:
        if p.exitcode != 0:
            is_success = False
    return is_success

def parser_commandline_args():
    desc = "NHSAP: A python pipeline  for NIPT non-human sequences analysis."

    cmdparse = argparse.ArgumentParser(description=desc)
    commands = cmdparse.add_subparsers(dest="command", title="NNSAP Commands")

    #For getting Non-human sequences
    gnhs_cmd = commands.add_parser('GNHS', help='Filter out Non-human sequence after alignmernt')
    gnhs_cmd.add_argument('-I', '--input', dest='input', metavar='fq', required=True, 
    					help='Raw fastq file containing human reads.')
    gnhs_cmd.add_argument('-b','--bwa',dest='bwa',metavar='bwa',
                        help='bwa path')
    gnhs_cmd.add_argument('-s','--samtools',dest='samtools',metavar='samtools',
                        help='samtools path')
    gnhs_cmd.add_argument('-bg','--bgzip',dest='bgzip',metavar='bgzip',
                        help='bgzip path')
    gnhs_cmd.add_argument('-tx','--tabix',dest='tabix',metavar='tabix',
                        help='tabix path')
    gnhs_cmd.add_argument('-si','--sample_id',dest='sample_id',metavar='sample_id',
                        help='sample_id')
    gnhs_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='Reference_fasta', required=True,
                        help='Input reference fasta file.')
    gnhs_cmd.add_argument('-c','--classification',dest='classification',metavar='classification',type=ast.literal_eval,
                        default=False,help='Read classification from bwa bam')
    gnhs_cmd.add_argument('-n','--nCPU',dest='nCPU',metavar='nCPU',
                        default='1',help='thread number for alignmernt')
    gnhs_cmd.add_argument('-O','--outdir',dest='outdir',metavar='fq',required=True,
    					help='Filtered non-human fq')

    #For Non-human sequences analysis
    gnha_cmd = commands.add_parser('GNHA', help='Non-human sequnce classification and profiling')
    gnha_cmd.add_argument('-I', '--input', dest='input', metavar='fq', required=True, 
    					help='non-human fq that has filered out human reads.')
    gnha_cmd.add_argument('-k2','--kraken2',dest='kraken2',metavar='kraken2',required=True,
    					help='kraken2 path')
    gnha_cmd.add_argument('-k2db','--kraken2database',dest='k2db',metavar='kraken2database',required=True,
    					help='kraken2 database path')
    	#kraken2 --db library fq.gz --quick --use-names --gzip-compressed --threads 1 --use-mpa-style --report  Kraken2.out --output kraken2.xls
    gnha_cmd.add_argument('-m2','--metaphlan2',dest='metaphlan2',metavar='metaphlan2',required=True,
    					help='metaphlan2 path')
    gnha_cmd.add_argument('-m2db','--metaphlan2database',dest='m2db',metavar='metaphlan2database',required=True,
    					help='metaphlan2 bowtie2 database path prefix')
    gnha_cmd.add_argument('-m2pkl','--metaphlan2pkl',dest='m2pkl',metavar='metaphlan2pkl database',required=True,
    					help='metaphlan2 pkl database path')
    	#metaphlan2.py fq.gz --input_type fastq --sample_id sampleid --bowtie2db ./db_v20/mpa_v20_m200 --bowtie2out bowtie2out.bz2 --mpa_pkl ./db_v20/mpa_v20_m200.pkl --nproc 1 --stat avg_g -t rel_ab_w_read_stats > MetaPhlAn2.out
    gnha_cmd.add_argument('-si','--sample_id',dest='sample_id',metavar='sample_id',
                        help='sample_id')
    gnha_cmd.add_argument('-O','--outdir',dest='outdir',metavar='outdir',required=True,
                        help='outdir of non-human analysis')
    gnha_cmd.add_argument('-p','--panel',dest='panel',metavar='panel',required=False,
                        help='non-human species panel')
    gnha_cmd.add_argument('-d','--db',dest='db',metavar='panel',required=False,
                        help='non-human references database')
    gnha_cmd.add_argument('-bs','--bamstats_summury',dest='bamstats_summury',metavar='panel',required=False,
                        help='mapping summary infomation to human reference genome')
    gnha_cmd.add_argument('-b','--bwa',dest='bwa',metavar='panel',required=False,
                        help='bwa path')
    gnha_cmd.add_argument('-n','--nCPU',dest='nCPU',metavar='nCPU',
                        default='1',help='thread number')
    gnha_cmd.add_argument('-s','--samtools',dest='samtools',metavar='panel',required=False,
                        help='samtools path')
    gnha_cmd.add_argument('-bg','--bgzip',dest='bgzip',metavar='panel',required=False,
                        help='bgzip path')
    gnha_cmd.add_argument('-tx','--tabix',dest='tabix',metavar='panel',required=False,
                        help='tabix path')

    return cmdparse.parse_args()

def main():
    start_time = time.time()
    runner = {
        'GNHS': gnhs,
        'GNHA': gnha,
    }

    args = parser_commandline_args()
    print(args)
    sys.stderr.write('\n** %s Start at %s **\n\n' % (args.command, time.asctime()))

    is_success = runner[args.command](args)

    elapsed_time = time.time() - start_time
    if is_success:
        sys.stderr.write('** %s done at %s, %d seconds elapsed **\n' % (
            args.command, time.asctime(), elapsed_time))
    else:
        sys.stderr.write('[ERROR] Catch some exception on %s, so "%s" is not done, %d seconds elapsed\n' % (
            time.asctime(), args.command, elapsed_time))
        sys.exit(1)

if __name__ == '__main__':
    main()
