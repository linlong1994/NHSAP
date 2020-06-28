
import os 
import sys
import glob
import re
from . import process_runner
from multiprocessing import Process

this_file = os.path.realpath(__file__)
this_file_dir = os.path.dirname(this_file)

#检查软件是否存在
def _precondition(software,softwarePath,warningString):

	if softwarePath == None:
		if software in os.environ:
			return True
		else:
			sys.stderr.write("[ERROR] There no %s in your environ variable please set --%s.\n\n" % (software,warningString))
			sys.exit(1)
	else:
		if os.path.exists(softwarePath):
			return True
		else:
			if software in os.environ:
				sys.stderr.write("[Warining] %s No such file or directory, and we have set %s from environment variable.\n\n" % (softwarePath,software))
				return True
			else:
				sys.stderr.write("[ERROR] %s No such file or directory,please re-set --%s.\n\n" % (softwarePath,warningString))
				sys.exit(1)

#检查index referecne是否存在
def _reference_check(software,reference_path,warningString):
	if reference_path == None or (not os.path.exists(reference_path)):
		sys.stderr.write("[ERROR] %s No such file or directory, please re-set --%s.\n\n " % (reference_path,warningString))
		sys.exit(1)
	else:
		if software == 'samtools':
			samtools_index = reference_path + '.fai'
			if not os.path.exists(samtools_index):
				sys.stderr.write("[ERROR] %s No samtools indexed reference file.\n\n" % (samtools_index))
				sys.exit(1)
			else:
				return True
		if software == 'bwa':
			bwa_index1 = reference_path + '.amb'
			bwa_index2 = reference_path + '.ann'
			bwa_index3 = reference_path + '.bwt'
			bwa_index4 = reference_path + '.pac'
			bwa_index5 = reference_path + '.sa'
			bwa_index_exists = os.path.exists(bwa_index1) and os.path.exists(bwa_index2) and os.path.exists(bwa_index3) and os.path.exists(bwa_index4) and os.path.exists(bwa_index5)
			if bwa_index_exists:
				return True
			else:
				sys.stderr.write("[ERROR] %s No bwa indexed reference file.\n\n" % (reference_path))
				sys.exit(1)
	return True

def _align(fq,sample_id,outdir,faidxed_reference,indexed_reference,bwa,samtools,bgzip,tabix,classification=False,stats=False,cpu=1):	
	sys.stderr.write(" %s aln -e 10 -t %s -i 5 -q 0 %s %s > %s/%s.sai\n" % (bwa,cpu,indexed_reference,fq,outdir,sample_id))
	os.system("time %s aln -e 10 -t %s -i 5 -q 0 %s %s > %s/%s.sai" % (bwa,cpu,indexed_reference,fq,outdir,sample_id))
	sys.stderr.write('%s samse -r "@RG\\tID:%s\\tPL:NA\\tSM:%s" %s %s/%s.sai %s | %s view -h -Sb - > %s/%s.bam \n' % (bwa,sample_id,sample_id,indexed_reference,outdir,sample_id,fq,samtools,outdir,sample_id))
	os.system('time %s samse -r "@RG\\tID:%s\\tPL:NA\\tSM:%s" %s %s/%s.sai %s | %s view -h -Sb - > %s/%s.bam && echo "** bwa done **" ' % (bwa,sample_id,sample_id,indexed_reference,outdir,sample_id,fq,samtools,outdir,sample_id))
	sys.stderr.write('time %s sort -@ 8 -O bam -o %s/%s.sorted.bam %s/%s.bam \n' % (samtools,outdir,sample_id,outdir,sample_id,) )
	os.system('time %s sort -@ 8 -O bam -o %s/%s.sorted.bam %s/%s.bam && echo "** bam sorted done **"' % (samtools,outdir,sample_id,outdir,sample_id,) )
	sys.stderr.write('time %s rmdup %s/%s.sorted.bam %s/%s.sorted.rmdup.bam \n' % (samtools,outdir,sample_id,outdir,sample_id) )
	os.system('time %s rmdup %s/%s.sorted.bam %s/%s.sorted.rmdup.bam && echo "** rmdup done **"' % (samtools,outdir,sample_id,outdir,sample_id) )
	sys.stderr.write('time %s index %s/%s.sorted.rmdup.bam \n' % (samtools,outdir,sample_id) )
	
	os.system('time %s index %s/%s.sorted.rmdup.bam && echo "** index done **"' % (samtools,outdir,sample_id) )
	sys.stderr.write('time %s mpileup --reference %s --min-MQ 0 --min-BQ 10 %s/%s.sorted.rmdup.bam | %s -f > %s/%s.mpileup.gz && %s -f -b 2 -e 2 %s/%s.mpileup.gz \n' % (samtools,faidxed_reference,outdir,sample_id,bgzip,outdir,sample_id,tabix,outdir,sample_id) )
	os.system('time %s mpileup --reference %s --min-MQ 0 --min-BQ 10 %s/%s.sorted.rmdup.bam | %s -f > %s/%s.mpileup.gz && %s -f -b 2 -e 2 %s/%s.mpileup.gz && echo "** bgzip mpileup done **"' % (samtools,faidxed_reference,outdir,sample_id,bgzip,outdir,sample_id,tabix,outdir,sample_id) )
	os.system('rm -rf %s/%s.sai %s/%s.bam %s/%s.sorted.bam ' % (outdir,sample_id,outdir,sample_id,outdir,sample_id) )
	sorted_rmdup_bam = outdir+'/'+sample_id+'.sorted.rmdup.bam'
	if classification:
		reads_classification_file = _read_classification_from_bwa_bam(samtools,sorted_rmdup_bam,sample_id,outdir)
		# 'sample_id.read_classification.txt'
	if stats:
		sys.stderr.write('%s stats %s > %s\n' % (samtools,sorted_rmdup_bam,sorted_rmdup_bam+'stats') )
		os.system('%s stats %s > %s' % (samtools,sorted_rmdup_bam,sorted_rmdup_bam+'stats') )
		bamstatScript = this_file_dir + '/external/bamstats.pl'
		os.system('perl %s %s %s > %s' % (bamstatScript,sorted_rmdup_bam+'stats',sample_id,sorted_rmdup_bam+'stats.summary.txt'))
	return sorted_rmdup_bam

def _align(fq,sample_id,outdir,faidxed_reference,indexed_reference,bwa,samtools,bgzip,tabix,classification=False,stats=False,cpu=1):	
	sys.stderr.write(" %s single -t %s %s %s > %s/%s.bam" % (bwa,cpu,indexed_reference,fq,outdir,sample_id))
	os.system("time %s single -t %s %s %s > %s/%s.bam" % (bwa,cpu,indexed_reference,fq,outdir,sample_id))
	sys.stderr.write('time %s sort -@ 8 -O bam -o %s/%s.sorted.bam %s/%s.bam \n' % (samtools,outdir,sample_id,outdir,sample_id,) )
	os.system('time %s sort -@ 8 -O bam -o %s/%s.sorted.bam %s/%s.bam && echo "** bam sorted done **"' % (samtools,outdir,sample_id,outdir,sample_id,) )
	sys.stderr.write('time %s rmdup %s/%s.sorted.bam %s/%s.sorted.rmdup.bam \n' % (samtools,outdir,sample_id,outdir,sample_id) )
	os.system('time %s rmdup %s/%s.sorted.bam %s/%s.sorted.rmdup.bam && echo "** rmdup done **"' % (samtools,outdir,sample_id,outdir,sample_id) )
	sys.stderr.write('time %s index %s/%s.sorted.rmdup.bam \n' % (samtools,outdir,sample_id) )
	
	os.system('time %s index %s/%s.sorted.rmdup.bam && echo "** index done **"' % (samtools,outdir,sample_id) )
	sys.stderr.write('time %s mpileup --reference %s --min-MQ 0 --min-BQ 10 %s/%s.sorted.rmdup.bam | %s -f > %s/%s.mpileup.gz && %s -f -b 2 -e 2 %s/%s.mpileup.gz \n' % (samtools,faidxed_reference,outdir,sample_id,bgzip,outdir,sample_id,tabix,outdir,sample_id) )
	os.system('time %s mpileup --reference %s --min-MQ 0 --min-BQ 10 %s/%s.sorted.rmdup.bam | %s -f > %s/%s.mpileup.gz && %s -f -b 2 -e 2 %s/%s.mpileup.gz && echo "** bgzip mpileup done **"' % (samtools,faidxed_reference,outdir,sample_id,bgzip,outdir,sample_id,tabix,outdir,sample_id) )
	os.system('rm -rf %s/%s.sai %s/%s.bam %s/%s.sorted.bam ' % (outdir,sample_id,outdir,sample_id,outdir,sample_id) )
	sorted_rmdup_bam = outdir+'/'+sample_id+'.sorted.rmdup.bam'
	if classification:
		reads_classification_file = _read_classification_from_bwa_bam(samtools,sorted_rmdup_bam,sample_id,outdir)
		# 'sample_id.read_classification.txt'
	if stats:
		sys.stderr.write('%s stats %s > %s\n' % (samtools,sorted_rmdup_bam,sorted_rmdup_bam+'stats') )
		os.system('%s stats %s > %s' % (samtools,sorted_rmdup_bam,sorted_rmdup_bam+'stats') )
		bamstatScript = this_file_dir + '/external/bamstats.pl'
		os.system('perl %s %s %s > %s' % (bamstatScript,sorted_rmdup_bam+'stats',sample_id,sorted_rmdup_bam+'stats.summary.txt'))
	return sorted_rmdup_bam

def _single_referecen_abundance(reads_classification_file,bamstats_summury,reference_fai,species_name):

	#'sample_id.abundance.txt'
	abundanceFile = reads_classification_file.replace('read_classification','abundance')
	getAbunceScript = this_file_dir + '/external/getAbudance.pl'

	sys.stderr.write('perl %s %s %s %s %s > %s\n' % (getAbunceScript,bamstats_summury,reads_classification_file,reference_fai,species_name,abundanceFile))
	os.system('perl %s %s %s %s %s > %s' % (getAbunceScript,bamstats_summury,reads_classification_file,reference_fai,species_name,abundanceFile))
	return abundanceFile

def _algin_multiReference(fq,sample_id,outdir,species_list,reference_list,bamstats_summury,bwa,samtools,bgzip,tabix):

	for i in range(len(reference_list)):
		faidxed_reference_precondition = _reference_check('samtools',reference_list[i],'reference')
		bwa_reference_precondition = _reference_check('bwa',reference_list[i],'reference')
		reference_fai = reference_list[i] + '.fai'
		species_name =  species_list[i]
		if faidxed_reference_precondition and bwa_reference_precondition:
			final_outdir = outdir + '/' + species_name
			if not os.path.exists(final_outdir):
				os.makedirs(final_outdir)
			sortBam = _align(fq,sample_id,final_outdir,reference_list[i],reference_list[i],bwa,samtools,bgzip,tabix,classification=True)
			reads_classification_file = sortBam.replace('sorted.rmdup.bam','read_classification.txt')
			abundanceFile = _single_referecen_abundance(reads_classification_file,bamstats_summury,reference_fai,species_name)
	return






def _unmapped_fq_from_bam(samtools,bam,sample_id,outdir):

	unmapped_fq_from_bam_perl_script = this_file_dir + '/external/unmapped_FQ_from_bam.pl'
	#unmapped_fq_from_bam($samtools,$bam,$sample_id,$outdir)
	os.system('time perl %s %s %s %s %s && echo "**Unmmapped reads grasped in FQ done **"' % (unmapped_fq_from_bam_perl_script,samtools,bam,sample_id,outdir) )

	return outdir+'/'+sample_id+'.fq.gz'

def _read_classification_from_bwa_bam(samtools,bam,sample_id,outdir):

	read_classification_from_bwa_bam_perl_script = this_file_dir + '/external/read_classification_from_bwa_bam.pl'
	#read_classification_from_bwa_bam($samtools,$bam,$sample_id,$outdir);
	os.system('time perl %s %s %s %s %s && echo "**Reads classification done **"' % (read_classification_from_bwa_bam_perl_script,samtools,bam,sample_id,outdir))
	reads_classification_file = outdir + '/' + sample_id + '.read_classification.txt'
	return reads_classification_file

class GnhsRunner(object):

	def __init__(self,args):

		self.inputfq = args.input
		self.bwa = args.bwa
		self.samtools = args.samtools
		self.reference = args.referencefile 
		self.outdir = args.outdir
		self.bgzip = args.bgzip
		self.tabix = args.tabix
		self.sample_id = args.sample_id
		self.nCPU = args.nCPU
		self.classification = args.classification

	def gnhs_runner(self):
		bwa_precondition = _precondition('bwa',self.bwa,'bwa')
		samtools_precondition = _precondition('samtools',self.samtools,'samtools')
		bgzip_precondition = _precondition('bgzip',self.bgzip,'bgzip')
		tabix_precondition = _precondition('tabix',self.tabix,'tabix')
		faidxed_reference_precondition = _reference_check('samtools',self.reference,'reference')
		bwa_reference_precondition = _reference_check('bwa',self.reference,'reference')

		if bwa_precondition and samtools_precondition and bgzip_precondition and tabix_precondition and faidxed_reference_precondition and bwa_reference_precondition:
			# sys.stderr.write(self.inputfq,self.sample_id,self.output,self.reference,self.reference,self.bwa,self.samtools,self.bgzip,self.tabix)
			if not os.path.exists(self.outdir):
				os.makedirs(self.outdir)
			sorted_rmdup_bam = _align(self.inputfq,self.sample_id,self.outdir,self.reference,self.reference,self.bwa,self.samtools,self.bgzip,self.tabix,self.classification,stats=True,cpu=self.nCPU)
			unmapped_fq = _unmapped_fq_from_bam(self.samtools,sorted_rmdup_bam,self.sample_id,self.outdir)
			# if self.classification:
			# 	_read_classification_from_bwa_bam(self.samtools,sorted_rmdup_bam,self.sample_id,self.output)
			return True


def _metaphlan2(metaphlan2,fq,sample_id,bowtie2db,mpa_pkl,outdir,nproc=1,stat='avg_g',t='rel_ab_w_read_stats',input_type='fastq'):
	bowtie2out = outdir + '/' + sample_id + '.bowtie2out.bz2'
	metaphlan2out = outdir + '/' + sample_id + '.MetaPhlAn2.out'
	metaphlan2outOrder = outdir + '/' + sample_id + '.MetaPhlAn2.out.order'

	sys.stderr.write('time %s %s --input_type %s --sample_id %s --bowtie2db %s --bowtie2out %s --mpa_pkl %s --nproc %s --stat %s -t %s > %s \n' % (metaphlan2,fq,input_type,sample_id,bowtie2db,bowtie2out,mpa_pkl,nproc,stat,t,metaphlan2out))
	os.system('time %s %s --input_type %s --sample_id %s --bowtie2db %s --bowtie2out %s --mpa_pkl %s --nproc %s --stat %s -t %s > %s && echo "**metaphlan2 done **" '% (metaphlan2,fq,input_type,sample_id,bowtie2db,bowtie2out,mpa_pkl,nproc,stat,t,metaphlan2out) )
	#what 
	extracrScript = this_file_dir + '/external/MetaPhlAn2/extract.taxo.addStrain.splitEukaryota.stat.pl' 
	refFungi = this_file_dir + '/external/MetaPhlAn2/fungi.cut2.xls'
	refParasite = this_file_dir + '/external/MetaPhlAn2/parasite.cut2.xls'
	sys.stderr.write('perl %s %s %s %s %s > %s\n' % (extracrScript,refFungi,refParasite,metaphlan2out,outdir,metaphlan2outOrder))
	os.system('perl %s %s %s %s %s > %s' % (extracrScript,refFungi,refParasite,metaphlan2out,outdir,metaphlan2outOrder) )
	return metaphlan2out

def _db_check(software,db_path,warningString):

	if software == 'bowtie':
		files = glob.glob(db_path + '*.bt2')
		if len(files) >= 1:
			return True
		else:
			sys.stderr.write("[ERROR] %s No %s database. please re-set --%s \n\n" % (db_path,software,warningString))
			sys.exit(1)

def _metaphlan2run(metaphlan2,fq,sample_id,bowtie2db,mpa_pkl,outdir,nproc=1,stat='avg_g',t='rel_ab_w_read_stats',input_type='fastq'):

	metaphlan2_precondition = _precondition('MetaPhlAn2',metaphlan2,'metaphlan2')
	fq_precondition = _reference_check('fastq',fq,'input')
	bowtie2db_precondition = _db_check('bowtie',bowtie2db,'metaphlan2database')
	mpa_pkl_precondition = _reference_check('mpa_pkl',mpa_pkl,'mpa_pkl')

	if metaphlan2_precondition and fq_precondition and bowtie2db_precondition and mpa_pkl_precondition:
		sys.stderr.write('metaphlan2run begins\n')
		metaphlan2out = _metaphlan2(metaphlan2,fq,sample_id,bowtie2db,mpa_pkl,outdir,nproc=1,stat='avg_g',t='rel_ab_w_read_stats',input_type='fastq')
		return metaphlan2out
	else:
		sys.stderr.write('[ERROR] Catch some exception on metaphlan2, so "metaphlan2" is not done \n')
		sys.exit(1)

def _kraken2(kraken2,fq,sample_id,kraken2db,outdir,thread=1,):
	kraken2report = outdir + '/' + sample_id + '.Kraken2.out'
	kraken2output = outdir + '/' + sample_id + '.Kraken2.xls'
	kraken2reportReAbun = outdir + '/' + sample_id + '.Kraken2.out.reAbun.txtt'
	kraken2SpeciesReAbun =  outdir + '/Species.reAbu.txt'
	kraken2SpeciesReAbunKingdom = outdir + '/Species.addKingdom.reAbu.txt'

	os.system('time %s --db %s %s --quick --use-names --gzip-compressed --threads %s --use-mpa-style --report %s --output %s ' % (kraken2,kraken2db,fq,thread,kraken2report,kraken2output) )
	statSampleScript = this_file_dir + '/external/Kraken2/Stat_Sample_Species_Kraken2.pl'
	os.system('perl %s -f %s -o %s -p %s' % (statSampleScript,kraken2report,outdir,sample_id) )
	extracrScript = this_file_dir + '/external/Kraken2/kraken2.extract.species.pl'
	os.system('perl %s %s >%s 2>%s' % (extracrScript,kraken2reportReAbun,kraken2SpeciesReAbun,kraken2SpeciesReAbunKingdom) )

	return kraken2report

def _kraken2run(kraken2,fq,sample_id,kraken2db,outdir,thread=1):
	kraken2_preconditon = _precondition('kraken2',kraken2,'kraken2')
	fq_precondition = _reference_check('fastq',fq,'input')
	kraken2db_precondition = _reference_check('kraken2',kraken2db,'kraken2')
	if kraken2_preconditon and kraken2db_precondition and fq_precondition:
		kraken2out = _kraken2(kraken2,fq,sample_id,kraken2db,outdir,thread)
		return kraken2out
	else:
		sys.stderr.write('[ERROR] Catch some exception on kraken2, so "kraken2" is not done \n')
		sys.exit(1)

def _get_species_from_metaphlan2out(metaphlan2out,relative_abundance_threshold=0,coverage_threshold=0,estimated_number_of_reads_threshold=0):
	metaphlan2Species = []
	if relative_abundance_threshold+coverage_threshold+estimated_number_of_reads_threshold > 0:
		for line in open(metaphlan2out):
			#clade_name     relative_abundance      coverage        average_genome_length_in_the_clade      estimated_number_of_reads_from_the_clade
			if 's__' in line: 
				line = line.strip('\n')
				col = line.split('\t')
				if (float(col[1]) > relative_abundance_threshold) and (float(col[2]) >coverage_threshold) and (float(col[4]) > estimated_number_of_reads_threshold):
					speciesName = re.search('s__([^|\t]*)\|?',col[0]).group(1)
					if 'unclassified' not in speciesName:
						metaphlan2Species.append(speciesName)
	else:
		for line in open(metaphlan2out):
			if 's__' in line:
				speciesName = re.search('s__([^|\t]*)\|?',line).group(1)
				if 'unclassified' not in speciesName:
					metaphlan2Species.append(speciesName)
	return metaphlan2Species

def _get_species_from_kraken2out(kraken2out,estimated_number_of_reads_threshold=0):
	kraken2Species = []
	if estimated_number_of_reads_threshold > 0:
		for line in open(kraken2out):
			#clade_name estimated_number_of_reads_threshold
			if 's__' in line: 
				line = line.strip('\n')
				col = line.split('\t')
				if float(col[1]) > estimated_number_of_reads_threshold:
					speciesName = re.search('s__([^|\t]*)\|?',col[0]).group(1)
					if 'unclassified' not in speciesName:
						speciesName = speciesName.replace(' ','_')
						speciesName = speciesName.replace('.','')
						speciesName = speciesName.replace('/','_')
						speciesName = speciesName.replace('-','_')
						speciesName = speciesName.replace("'",'')
						speciesName = speciesName.replace("#",'')
						speciesName = speciesName.replace('(','')
						speciesName = speciesName.replace(')','')
						speciesName = speciesName.replace(':','')
						kraken2Species.append(speciesName)
	else:
		for line in open(kraken2out):
			if 's__' in line:
				speciesName = re.search('s__([^|\t]*)\|?',line).group(1)
				if 'unclassified' not in speciesName:
					speciesName = speciesName.replace(' ','_')
					speciesName = speciesName.replace('.','')
					speciesName = speciesName.replace('/','_')
					speciesName = speciesName.replace('-','_')
					speciesName = speciesName.replace("'",'')
					speciesName = speciesName.replace("#",'')
					speciesName = speciesName.replace('(','')
					speciesName = speciesName.replace(')','')
					speciesName = speciesName.replace(':','')
					kraken2Species.append(speciesName)
	return kraken2Species

def _merge(metaphlan2Species,kraken2Species,db,panel):
	#panel : species names 
	#db : organism,refseq_category,assembly_level,ftp_genome
	panel_species = []
	for line in open(panel):
		line = line.strip('\n')
		col = line.split('\t')
		panel_species.append(col[0])
	db_reference = {}
	for line in open(db):
		line = line.strip('\n')
		col = line.split('\t')
		db_reference[col[0]] = col[-1]
	final_reference = []
	final_species = []
	speciesList = list(set(metaphlan2Species) | set(kraken2Species) | set(panel_species))
	for species in speciesList:
		if species in db_reference:
			final_reference.append(db_reference[species])
			final_species.append(species)
		else:
			sys.stderr.write("[Warining] No reference file for species: %s\n" % (species))
	return final_species,final_reference

def _div_list(ls,n):
	if not isinstance(ls,list) or not isinstance(n,int):
		return []
	ls_len = len(ls)
	if n<=0 or 0==ls_len:
		return []
	if n > ls_len:
		return []
	elif n == ls_len:
		return [[i] for i in ls]
	else:
		j = ls_len//n
		k = ls_len%n
		### len(ls) = (n-1)*j + (j+k)
		ls_return = []
		for i in range(0,(n-1)*j,j):
			ls_return.append(ls[i:i+j])
		ls_return.append(ls[(n-1)*j:])
		return ls_return

def _cat_files(file_list,output):
	O = open(output,'wt')
	for file in file_list:
		t = open(file)
		for line in t:
			O.write(line)
		t.close()
	O.close()
	return 

class GnhaRunner(object):

	def __init__(self,args):
		self.metaphlan2 = args.metaphlan2
		self.m2db = args.m2db
		self.m2pkl = args.m2pkl
		self.kraken2 = args.kraken2
		self.kraken2db = args.k2db
		self.fq = args.input
		self.sample_id = args.sample_id
		self.outdir = args.outdir
		self.panel = args.panel
		self.db = args.db
		self.bamstats_summury = args.bamstats_summury
		self.bwa = args.bwa
		self.samtools = args.samtools
		self.bgzip = args.bgzip
		self.tabix = args.tabix
		self.nCPU = args.nCPU

	def gnha_runner(self):

		metaphlan2outDir = self.outdir + '/MetaPhlan2'
		if not os.path.exists(metaphlan2outDir):
			os.makedirs(metaphlan2outDir)
		kraken2outDir = self.outdir + '/Kraken2'
		if not os.path.exists(kraken2outDir):
			os.makedirs(kraken2outDir)
		nhsapDir = self.outdir + '/nhsap'
		if not os.path.exists(nhsapDir):
			os.makedirs(nhsapDir)
		final_output = self.outdir + '/' + self.sample_id + '.abundance.txt'
		sys.stderr.write('metaphlan2 begins\n')
		#_metaphlan2run(metaphlan2,fq,sample_id,bowtie2db,mpa_pkl,outdir,nproc=1,stat='avg_g',t='rel_ab_w_read_stats',input_type='fastq')
		metaphlan2out = _metaphlan2run(self.metaphlan2,self.fq,self.sample_id,self.m2db,self.m2pkl,metaphlan2outDir,nproc=self.nCPU)
		sys.stderr.write('kraken2 begins\n')
		# _kraken2run(kraken2,fq,sample_id,kraken2db,outdir,thread=1)
		kraken2out = _kraken2run(self.kraken2,self.fq,self.sample_id,self.kraken2db,kraken2outDir,thread=self.nCPU)
		#_get_species_from_metaphlan2out(metaphlan2out,relative_abundance_threshold=0,coverage_threshold=0,estimated_number_of_reads_threshold=0)
		metaphlan2Species = _get_species_from_metaphlan2out(metaphlan2out,relative_abundance_threshold=0,coverage_threshold=0.00001,estimated_number_of_reads_threshold=10)
		#_get_species_from_kraken2out(kraken2out,estimated_number_of_reads_threshold=0)
		kraken2Species = _get_species_from_kraken2out(kraken2out,estimated_number_of_reads_threshold=0)

		final_species,final_reference = _merge(metaphlan2Species,kraken2Species,self.db,self.panel)
		print(final_species)

		final_reference_batch = _div_list(final_reference,int(self.nCPU))
		final_species_batch = _div_list(final_species,int(self.nCPU))

		processes = []
		abundance_batch_files = []
		for i in range(int(self.nCPU)):
			sub_dir = nhsapDir + '/batch_%s' % i
			if not os.path.exists(sub_dir):
				os.makedirs(sub_dir)
			abundance_batch_files += ['%s/%s/%s.abundance.txt' % (sub_dir,final_species_batch[i][u],self.sample_id) for u in range(len(final_species_batch[i])) ]
			#_algin_multiReference(fq,sample_id,outdir,species_list,reference_list,bamstats_summury,bwa,samtools,bgzip,tabix)
			processes.append(Process(target = _algin_multiReference,args=(
											self.fq,
											self.sample_id,
											sub_dir,
											final_species_batch[i],
											final_reference_batch[i],
											self.bamstats_summury,
											self.bwa,
											self.samtools,
											self.bgzip,
											self.tabix
											)))
		process_runner(processes)

		_cat_files(abundance_batch_files,final_output)

		return processes

