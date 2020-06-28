
import os 
import sys 
import math

def log100(x):
	import math
	return -3*(math.log(x * 1000000)  / math.log(0.000001))-3


def smartDownload(filelist,outdir,bwa=None,samtools=None):
	
	downloadFileList = []

	fileNum = len(filelist)
	dir_layer_num = math.ceil(log100(fileNum)) - 1
	sub_dirs = ['1'] * dir_layer_num
	fileCount = 0
	for file in filelist:
		filename = file.split('/')[-1]
		fileCount += 1
		if fileCount == 101:
			sub_dirs[-1] =str(int(sub_dirs[-1])+1)
		for x in range(len(sub_dirs)-1,-1,-1):
			if sub_dirs[x] == '101':
				sub_dirs[x] = '1'
				sub_dirs[x-1] = str(int(sub_dirs[x-1]) + 1)
		if sub_dirs:
			final_dir = outdir + '/' +  '/'.join(sub_dirs)
			if not os.path.exists(final_dir):
				os.makedirs(final_dir)
			final_file = final_dir + '/' + file_name
		#wget(file,final_file)
			if bwa:
				index(final_file,samtools,bwa)
			downloadFileList.append(final_file)
	return downloadFileList

def index(referecne,samtools,bwa):
	sys.stderr.write('%s index %s' % (bwa,referecne))
	os.system('%s index %s' % (bwa,referecne))
	sys.stderr.write('%s faidx %s' % (samtools,referecne))
	os.system('%s faidx %s' % (samtools,referecne))


def referenceDownload(assembly_report,bwa=None,samtools=None,filetype='fna',outdir='./'):

	if assembly_report:
		pass
	else:
		assembly_report = outdir + '/' + 'assembly_summary_refseq.txt'
		os.system('wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt %s' % assembly_report)

	summary = open(assembly_report)
	refseq_summary = {}
	final_organism = []
	final_ftp_genomes = []
	for line in summary:
		if line.startswith('#'):
			continue
		line = line.strip('\n')
		col = line.split('\t')
		organism_name = col[7]
		assembly_level = col[11] #Chromosome	Complete Genome	Contig	Scaffold
		refseq_category = col[4] #na	reference genome	representative genome
		ftp_genome = col[19]
		if organism_name not in refseq_summary:
			refseq_summary[organism_name] = {}
		if refseq_category not in refseq_summary[organism_name]:
			refseq_summary[organism_name][refseq_category] = {}
		if assembly_level not in refseq_summary[organism_name][refseq_category]:
			refseq_summary[organism_name][refseq_category][assembly_level] = []
		refseq_summary[organism_name][refseq_category][assembly_level].append(ftp_genome)	

	for organism in refseq_summary:
		if  'reference genome' in refseq_summary[organism]:
			refseq_category = 'reference genome'
			if 'Complete Genome' in refseq_summary[organism]['reference genome']:
				assembly_level = 'Complete Genome'
				ftp_genome = refseq_summary[organism]['reference genome']['Complete Genome'][0]
			elif 'Chromosome' in refseq_summary[organism]['reference genome']:
				assembly_level = 'Chromosome'
				ftp_genome = refseq_summary[organism]['reference genome']['Chromosome'][0]
			elif 'Scaffold' in refseq_summary[organism]['reference genome']:
				assembly_level = 'Scaffold'
				ftp_genome = refseq_summary[organism]['reference genome']['Scaffold'][0]
			elif 'Contig' in refseq_summary[organism]['reference genome']:
				assembly_level = 'Contig'
				ftp_genome = refseq_summary[organism]['reference genome']['Contig'][0]
		elif 'representative genome' in refseq_summary[organism]:
			refseq_category = 'representative genome'
			if 'Complete Genome' in refseq_summary[organism]['representative genome']:
				assembly_level = 'Complete Genome'
				ftp_genome = refseq_summary[organism]['representative genome']['Complete Genome'][0]
			elif 'Chromosome' in refseq_summary[organism]['representative genome']:
				assembly_level = 'Chromosome'
				ftp_genome = refseq_summary[organism]['representative genome']['Chromosome'][0]
			elif 'Scaffold' in refseq_summary[organism]['representative genome']:
				assembly_level = 'Scaffold'
				ftp_genome = refseq_summary[organism]['representative genome']['Scaffold'][0]
			elif 'Contig' in refseq_summary[organism]['representative genome']:
				assembly_level = 'Contig'
				ftp_genome = refseq_summary[organism]['representative genome']['Contig'][0]
		else:
			refseq_category = 'na'
			if 'Complete Genome' in refseq_summary[organism]['na']:
				assembly_level = 'Complete Genome'
				ftp_genome = refseq_summary[organism]['na']['Complete Genome'][0]
			elif 'Chromosome' in refseq_summary[organism]['na']:
				assembly_level = 'Chromosome'
				ftp_genome = refseq_summary[organism]['na']['Chromosome'][0]
			elif 'Scaffold' in refseq_summary[organism]['na']:
				assembly_level = 'Scaffold'
				ftp_genome = refseq_summary[organism]['na']['Scaffold'][0]
			elif 'Contig' in refseq_summary[organism]['na']:
				assembly_level = 'Contig'
				ftp_genome = refseq_summary[organism]['na']['Contig'][0]
		if filetype == 'fna':
			ftp_genome = ftp_genome + '/' + ftp_genome.split('/')[-1] + '_genomic.fna.gz'

		organism = organism.replace(' ','_')
		
		final_organism.append('%s\t%s\t%s' % (organism,refseq_category,assembly_level))
		final_ftp_genomes.append(ftp_genome)
	for index in range(len(final_organism)):
		print(final_organism[index],'\t',final_ftp_genomes[index])
	# downloadFiles = smartDownload(final_ftp_genomes,bwa,samtools)
	# for index in range(len(downloadFiles)):
	# 	sys.stderr.write('%s\t%s\n' % (final_organism[index],downloadFiles[index]))
		#sys.stderr.write('%s\t%s\t%s\t%s' % (organism,refseq_category,assembly_level,ftp_genome))


		# for catory in refseq_summary[organism]:
		# 	for level in refseq_summary[organism][catory]:
		# 		if len(refseq_summary[organism][catory][level]) > 1 :
		# 			sys.stderr.write('%s\t%s\t%s' % (organism,catory,level))
		# 			sys.stderr.write(refseq_summary[organism][catory][level])


def mannual():
	sys.stderr.write('''
python3 downloadReference.py 输入对应软件的database referencePanel 然后根据assembly_refseq下载所需文件 
	需要生成MetaPhlAn2/fungi.cut2.xls MetaPhlAn2/parasite.cut2.xls
	''')


if __name__ == '__main__':
	referenceDownload('./assembly_summary_refseq.txt')

