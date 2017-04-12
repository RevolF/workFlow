20160705
get to know tmap params well with real datasets
try to automize optimization using python utilizing CV
read paper about that optimized pipeline and comparation between different softwares
most important,try break point detections

using MLPA for detection of indels on chromosome
this method only confirms which intron has duplicated or deleted but without power to pridict which point
array comparative genome hybridization(array-CGH), scanning copy number variation in human genome, a LS algorithsm was introduced but not pratical

SVMerge from sanger provides better tools for detecting of breakpoints
	http://svmerge.sourceforge.net/#download
BreakDancer:
	https://sourceforge.net/projects/breakdancer/?source=recommended
these two may be better worked for next generation data
a cancer genetics reviews talks about all available informatics methods
	http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4441822/pdf/nihms-542544.pdf

several indel detection method
Pindel
GATK
VarScan2
Dindel: currently only for illumina data
Stanpy: illumina
SAMtools
breakdance

DMD/BMD project
using existing data from DMD/BMD clinical data at:
	/media/disk2/ywang/rawdata/rawdata_35/17panel/Auto_sn247560054_SN2-353-P35_Hi-Q-pooling-P17-33ZA1-20160630_546_824/basecaller_results
other important files are listed as:
	IDP17_bed
	/media/disk2/ybliu/project/monogenic/panel_design/bed2gene/IDP17_bed/0770341_Covered.bed

hg19.fasta
	/media/disk2/ybliu/project/WGS/data/hg19/ucsc.hg19.fasta

refGene.txt
	/home/ybliu/data/refGene.txt
	
	DMD:NM_004006
	
	ionadmin@192.168.0.35
	ionadmin
	
20160706
Pindel:
	Pindel can detect breakpoints of large deletions, medium sized insertions, inversions, tandem duplications and other structural variants at single-based resolution from next-gen sequence data.
	parameters explaination presented here:
		http://gmt.genome.wustl.edu/packages/pindel/quick-start.html
		http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
	for inserted sequence length in pindel config file: 250 were denoted to paired end sequencing method while for ion, this maybe smaller
	
Improving Indel Detection Specificity of the Ion Torrent PGM Benchtop Sequencer
	http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0045798
an example of getting indels from ion data in paper:
	Evaluation and optimisation of indel detection workflows for ion torrent sequencing of the BRCA1 and BRCA2 genes

variants caller on proton
	https://www.edgebio.com/variant-calling-ion-torrent-data

VCF file format:
	http://vcftools.sourceforge.net/specs.html

	
	
20160707-20160708
GATK tools tutorials can be found at:
	https://www.broadinstitute.org/gatk/guide/topic?name=tutorials
	documentations at:
	https://www.broadinstitute.org/gatk/guide/tooldocs/
	
multiple tools can be found at pengyuan\'s bin:
	47:/home/pyzhu/bin
	
	
commends:
	samtools view *.bam | less -S
	samtools index *.bam
	nohup samtools index IonXpress_001_rawlib.basecaller.bam >/dev/null 2>/dev/null &
	nohup pindel -f ../hg19_reference/ucsc.hg19.fasta -i bam_sample_config.txt -c chrX -T 20 -o 0707pindelOutput/ >nohupPindel.out &
	samtools mpileup -d 10000 -L 1000 -Q 7 -h 50 -o 10 -e 17 -m 4 -f ../hg19_reference/ucsc.hg19.fasta -g IonXpress_001_rawlib.basecaller.bam 
	samtools sort file.bam outputPrefix
	tmap index -vf ucsc.hg19.fasta
	tmap mapall -n 24 -f ucsc.hg19.fasta -r IonXpress_001_rawlib.basecaller.bam -v -Y -u -a 3 -s results001.bam -o 0 stage1 map4
	samtools sort -T ./samtmp/001sorted -o results001Sorted.bam results001.bam
	
	this following workflow should be a standard:
		nohup tmap mapall -f ucsc.hg19.fasta -r IonXpress_001_rawlib.basecaller.bam -s result001.bam stage1 map4 &
		nohup samtools sort -T ./tmp/re1 -o results_sorted.bam result001.bam &
		nohup samtools mpileup -f ucsc.hg19.fasta results_sorted.bam > results1.mpileup &
		nohup samtools mpileup -f ucsc.hg19.fasta -g results_sorted.bam -o results1mpileup.bcf &
			add a -g to ensure that the output is a bcf file forusing  bcftools
	
	
	
	
	
infos:
	http://bioinformatics.lofter.com/post/bffd5_a87dee
	http://blog.sciencenet.cn/blog-252888-731246.html
	for tmap params:
		http://129.130.90.13/ion-docs/TMAP-Modules.html
	
tips:
	for raw bam files from ion; first call tmap/bwa for mapping, then samtools for sorting, samtools then bcf
	a mpileup shell file could be found at: /media/disk2/ybliu/project/monogenic/tmp/WES_somatic_hlliu/WES_somatic_compare_160630
	use default params amap
	keep commends in a shell


0711
	creating java on .bashrc:
		export PATH=/home/ljzhang/bin/jdk1.8.0_91/bin:$PATH
		export PATH=/home/ljzhang/bin/jdk1.8.0_91/jre/bin:$PATH
		export JAVA_HOME=/home/ljzhang/bin/jdk1.8.0_91:$JAVA_HOME
		export CLASSPATH=/home/ljzhang/bin/jdk1.8.0_91/lib/tools.jar:$CLASSPATH
		export CLASSPATH=/home/ljzhang/bin/jdk1.8.0_91/lib/dt.jar:$CLASSPATH
	
	TVC parameter file: http://129.130.90.13/ion-docs/Example-Torrent-Variant-Caller-Parameter-File_51413889.html
	TVC commend line: http://129.130.90.13/ion-docs/The-Command-Line-Torrent-Variant-Caller.html
	
	fuzzy installing https://github.com/domibel/IonTorrent-VariantCaller
		write all cmds in a shell and excute
		this procedure requires boost C++ upport, which could be download and compiled from source-forge
		
	picard: https://broadinstitute.github.io/picard/
		https://github.com/broadinstitute/picard
		
	tvc usage:
		https://ioncommunity.thermofisher.com/docs/DOC-9504
	/home/ljzhang/local/ion_pvc_install.sh
	/home/ljzhang/sftwdl/boost_1_61_0.tar.gz
		./bootstrap || ./b2 
		http://blog.csdn.net/lixiang987654321/article/details/49894987
	this may require sudo for boost installing

0712
	change export method, solving ant and JAVA_HOME problems
		PATH="/home/ljzhang/bin/jdk1.8.0_91/bin${PATH:+:${PATH}}"; export PATH
	grep 'INDEL' /home/ljzhang/project/dmd/0707_bam_test/pipetest08/chrXbcf.out | less -S
	
	installing htsjdk for picard:
		https://github.com/samtools/htsjdk
			* What went wrong:
			A problem occurred configuring root project 'htsjdk'.
			> Could not resolve all dependencies for configuration ':classpath'.
			   > Could not download shadow.jar (com.github.jengelman.gradle.plugins:shadow:1.2.3)
				  > Could not get resource 'https://plugins.gradle.org/m2/com/github/jengelman/gradle/plugins/shadow/1.2.3/shadow-1.2.3.jar'.
					 > Could not GET 'https://plugins.gradle.org/m2/com/github/jengelman/gradle/plugins/shadow/1.2.3/shadow-1.2.3.jar'.
						> Remote host closed connection during handshake

	picard compilation:
		https://github.com/broadinstitute/picard/releases/tag/2.5.0
		https://github.com/broadinstitute/picard
		downloaded, always better to ask for new pathways
		/home/ljzhang/bin/picard-tools-2.5.0
		
	GATK runtime error is caused by java version
	
	in perl function, path is determined by the excution environment
	
	an effective way to kill all processes of a user
		killall -u ttlsa
		pgrep -u ttlsa | sudo xargs kill -9
		ps -ef | grep ttlsa | awk '{ print $2 }' | sudo xargs kill -9
		pkill -u ttlsa
	
	a bcf file specs is on 45 server at:
		/home/ljzhang/project/papaers/hts-specs
		
	use .pipelog to hold process infos to avoid repeats of processes
	future adding parameters using R
	currently focusing on the output of results, then come back for python subprocess management
	
	introduction to vcf 4.2 
		https://faculty.washington.edu/browning/intro-to-vcf.html
		

0713-0714
		java -jar VarScan somatic 
	use fastqc to check sequencing quality
	refer to a new model of finding breakpoints on genome
	GATK is a pipeline with algorithsms to do different work, which should be better regarded as a framework than a single software
		
	
	still use tmap for mapping, use bam as output, while applying picard for sort and mark duplications
		
	seems that sorting with samtools is faster than using picard
			
										##################################################
						############################### pipeline starts here #############################
										##################################################
		samtools index IonXpress_001_rawlib.basecaller.bam
		
		tmap mapall 
		-f ucsc.hg19.fasta 
		-r IonXpress_001_rawlib.basecaller.bam 
		-s tmap_mapall_IonXpress_001_rawlib.basecaller.bam 
		stage1 map4
		
		################################################################################################################
		# this following is used for bam sam converting
		samtools view -h 1.bam -o 1.sam			# adding -h for including headers in final sam output
		
		java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/picard-tools-2.5.0/picard.jar 
		SamFormatConverter 
		I=results_sorted_IonXpress_001_rawlib.basecaller.sam 
		O=picard_converted_001.bam
		
		################################################################################################################
		
		# picard sort
		java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/picard-tools-2.5.0/picard.jar 
		SortSam 
		I=tmap_mapall_IonXpress_001_rawlib.basecaller.bam 
		O=picard_sorted_IonXpress_001_rawlib.basecaller.bam 
		SO=coordinate
		#picard mark duplicates
		java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/picard-tools-2.5.0/picard.jar
		MarkDuplicates 
		I=picard_sorted_IonXpress_001_rawlib.basecaller.bam 
		O=picard_sorted_IonXpress_001_rawlib.basecaller_dedup.bam 
		M=picard_sorted_IonXpress_001_rawlib.basecaller_dedup.metrics 
		REMOVE_DUPLICATES=false 
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000
		
		##########################################
		#samtools depth, optional
		samtools depth picard_sorted_IonXpress_001_rawlib.basecaller_dedup.bam > qc_IonXpress_001_rawlib.basecaller_dedup.depth
		perl creat_dep_index.pl qc_IonXpress_001_rawlib.basecaller_dedup.depth
		perl depth_stat.pl qc_IonXpress_001_rawlib.basecaller_dedup.depth qc_IonXpress_001_rawlib.basecaller_dedup.depth.stat
		#
		
		# fastqc which is optional
		perl ~/bin/FastQC/fastqc -q -o ./fastqctmp picard_sorted_IonXpress_001_rawlib.basecaller_dedup.bam
		##########################################
		
		
		# GATK RealignerTargetCreator
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T RealignerTargetCreator 
		-R ucsc.hg19.fasta 
		-I picard_sorted_IonXpress_001_rawlib.basecaller_dedup.bam 
		-o picard_sorted_IonXpress_001_rawlib.basecaller_dedup.intervals 
		-nt 8 
		-known ./hg19/1000G_phase1.indels.hg19.sites.vcf 
		-known ./hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
		
		# GATK IndelRealigner
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T IndelRealigner 
		-R ucsc.hg19.fasta 
		-targetIntervals picard_sorted_IonXpress_001_rawlib.basecaller_dedup.intervals 
		-I picard_sorted_IonXpress_001_rawlib.basecaller_dedup.bam 
		-o picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.bam 
		-known ./hg19/1000G_phase1.indels.hg19.sites.vcf 
		-known ./hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
		
		# breakdancer
		
		#
		
		# GATK BaseRecalibrator
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T BaseRecalibrator 
		-R ucsc.hg19.fasta 
		-I picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.bam 
		-O picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.grp 
		-nct 8 
		-knownSites ./hg19/dbsnp_138.hg19.vcf 
		-knownSites ./hg19/1000G_phase1.indels.hg19.sites.vcf 
		-knownSites ./hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf 
		-bqsrBAQGOP 30
		
		# GATK PrintReads
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T PrintReads 
		-R ucsc.hg19.fasta 
		-I icard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.bam 
		-BQSR picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.grp 
		-o picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.recal.bam 
		-nct 8
		
		# GATK HaplotypeCaller
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T HaplotypeCaller 
		-R ucsc.hg19.fasta 
		-I picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.recal.bam 
		-o picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.recal.raw.vcf 
		-nct 4 
		-D ./hg19/dbsnp_138.hg19.vcf 
		-stand_call_conf 30 
		-stand_emit_conf 10 
		-rf BadCigar
		
		# GATK VariantRecalibrator SNP
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T VariantRecalibrator 
		-R ucsc.hg19.fasta 
		-input picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.recal.raw.vcf 
		-recalFile picard_sorted_IonXpress_001_rawlib.basecaller.snp.recal 
		-tranchesFile picard_sorted_IonXpress_001_rawlib.basecaller.snp.tranches 
		-rscriptFile picard_sorted_IonXpress_001_rawlib.basecaller.snp.plot.R 
		-mode SNP 
		-nt 8 
		-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP 
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./hg19/hapmap_3.3.hg19.sites.vcf 
		-resource:omni,known=false,training=true,truth=true,prior=12.0 ./hg19/1000G_omni2.5.hg19.sites.vcf 
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 ./hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf 
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./hg19/dbsnp_138.hg19.vcf 
		
		# GATK ApplyRecalibration SNP
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T ApplyRecalibration 
		-R ucsc.hg19.fasta 
		-input picard_sorted_IonXpress_001_rawlib.basecaller_dedup.realn.recal.raw.vcf 
		-recalFile picard_sorted_IonXpress_001_rawlib.basecaller.snp.recal 
		-tranchesFile picard_sorted_IonXpress_001_rawlib.basecaller.snp.tranches 
		-o picard_sorted_IonXpress_001_rawlib.basecaller.snp.filtered.vcf 
		-nt 8 
		-ts_filter_level 99.5 
		-mode SNP
		
		# GATK VariantRecalibrator INDEL
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T VariantRecalibrator 
		-R ucsc.hg19.fasta 
		-input picard_sorted_IonXpress_001_rawlib.basecaller.snp.filtered.vcf 
		-recalFile picard_sorted_IonXpress_001_rawlib.basecaller.indel.recal 
		-tranchesFile picard_sorted_IonXpress_001_rawlib.basecaller.indel.tranches 
		-rscriptFile picard_sorted_IonXpress_001_rawlib.basecaller.indel.plot.R 
		-mode INDEL 
		-nt 8 
		-resource:mills,known=true,training=true,truth=true,prior=12.0 ./hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf 
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./hg19/dbsnp_138.hg19.vcf 
		--maxGaussians 4 
		-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum 
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
		
		# GATK ApplyRecalibration INDEL
		~/usr/bin/java -Xmx40g -Djava.io.tmpdir=./javatmp/ -jar ~/bin/GATK/GenomeAnalysisTK/GenomeAnalysisTK.jar 
		-T ApplyRecalibration 
		-R ucsc.hg19.fasta 
		-input picard_sorted_IonXpress_001_rawlib.basecaller.snp.filtered.vcf 
		-recalFile picard_sorted_IonXpress_001_rawlib.basecaller.indel.recal 
		-tranchesFile picard_sorted_IonXpress_001_rawlib.basecaller.indel.tranches 
		-o picard_sorted_IonXpress_001_rawlib.basecaller.snp.indel.filtered.vcf 
		-nt 8 
		-mode INDEL 
		-ts_filter_level 99.0
		
		
										##################################################
						############################### pipeline ends here ###############################
										##################################################
	
		
		
		#### samtools calling for read depth and specific regions
		
		tmap mapall -f ucsc.hg19.fasta -r IonXpress_001_rawlib.basecaller.bam -s result001.bam stage1 map4
		samtools sort -T ./tmp/re1 -o results_sorted.bam result001.bam
		samtools depth results_sorted_IonXpress_001_rawlib.basecaller.bam -r chr1:1600000-1700000
		# using samtools for target genome intervals
		samtools view results_sorted_IonXpress_001_rawlib.basecaller.bam chr1:1600000-1700000
		
		####
	
	
	
	
	set JAVA_HOME to root dir of jdk rather than jre to avoid path conflict
	installing maven, use mvn install for pom.xml
	for installing GATK.3.6.tar.bz2, use unpack cmd of tar xjf GenomeAnalysisTK-3.6.tar.bz2, this would produce a gatk.jar file rather than the unpacking file of .jar (which is also a packed file)
		https://www.broadinstitute.org/gatk/guide/article?id=2899
	
	stick to official guide
	
	breakdancer:
		http://breakdancer.sourceforge.net/
		https://github.com/genome/breakdancer/blob/master/INSTALL.md
		now is compiled and included in PATH using cmd "breakdancer-max" to invoke
		
	hg19 sv files used for GATK could be found at:
		ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
	
	
0718
	samtools manual page:
		http://www.htslib.org/doc/samtools.html
		http://www.bbioo.com/lifesciences/40-112919-1.html
		
	be careful of space ident
	tvc is available on 47 server: /home/pyzhu/bin/tvc
	
	for gatk_pipeline.py, pipes proceeds to gatk haplotype
	
0719
	different file format:
		https://genome.ucsc.edu/FAQ/FAQformat.html
	using tvc pipeline on ionadmin@192.168.0.35
		/results/plugins/variantCaller/variant_caller_pipeline.py -i $out_dir/$sample_name/$sample_name/.fixed.bam -o $out_dir/$sample_name/$variant -B /results/plugins/variantCaller/ -n 10 -r $ref -b $target -p /results/plugins/KD_plugin/hotspot/germline_low_stringency_protons2.json -s /results/plugins/KD_plugin/hotspot/kd.hotspot.vcf
	remember to build an index for target bam before passing to tvc
		python /results/plugins/variantCaller/variant_caller_pipeline.py -i IonXpress_001_rawlib.basecaller.sorted.bam -o ./test/1 -B /results/plugins/variantCaller/ -n 10 -r ~/ljzhang/data/hg19/ucsc.hg19.fasta -p /results/plugins/KD_plugin/hotspot/germline_low_stringency_proton2.json -s /results/plugins/KD_plugin/hotspot/kd.hotspot.vcf
	
	
	
0720
	prepare for vcf info analysis.
	prepare for bam info analysis
	
	python /home/ljzhang/bin/variantCaller/variant_caller_pipeline.py -i /home/ljzhang/project/dmd/0707_bam_test/test_gatk_pipe_TMP/IonXpress_001_rawlib.basecaller.sorted.bam -o /home/ljzhang/project/dmd/0707_bam_test/test_gatk_pipe_TMP/test_tvc/IonXpress_001_rawlib.basecall -B /home/ljzhang/bin/variantCaller/ -n 4 -r /home/ljzhang/data/hg19/ucsc.hg19.fasta -p /home/ljzhang/data/KD_plugin/hotspot/germline_low_stringency_proton2.json -s /home/ljzhang/data/KD_plugin/hotspot/kd.hotspot.vcf
	
	tvc pipeline finished
	
	learn how to analysis a processed vcf file
	
	CIGAR and FLAG and QUALITY
	
	use CIGAR for extending sequence and finding breakpoint pannel
	a precise breakpoint would not be necessary, cutting human reference genome to 1000bp pannel from site to site and try to get an enrichment curve
	use bam file as raw input, first loose parameters of filter, by letting more reads in for tunning.
	for DMD project, we do not need to consider inversion, only get results out ASAP
	
	
0721
	for CIGAR string explaination:
		http://davetang.org/wiki/tiki-index.php?page=SAM
	
		Operator	Description
		D	Deletion; the nucleotide is present in the reference but not in the read
		H	Hard Clipping; the clipped nucleotides are not present in the read.
		I	Insertion; the nucleotide is present in the read  but not in the rference.
		M	Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
		N	Skipped region; a region of nucleotides is not present in the read
		P	Padding; padded area in the read and not in the reference
		S	Soft Clipping;  the clipped nucleotides are present in the read
		X	Read Mismatch; the nucleotide is present in the reference
		=	Read Match; the nucleotide is present in the reference
		
	
	using perl for parsing sam/bam:
		http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm
	biopython:
		http://biopython.org/DIST/docs/tutorial/Tutorial.html
	
	tmap->samtools sort -> samtools index -> samtools cut 
		knowgene downloaded at: /home/ljzhang/data/knownGene.txt
		samtools view IonXpress_001_rawlib.basecaller.sorted.bam chrX > IonXpress_001_rawlib.basecaller.chrX.sorted.sam
		
	python ./tmpscript/gethg19Chr.py
		output chrX to one single string
		using soft padding as main source of information
	
	DMD:NM_004006 transcript
	
	using python for getting desired line number:
	avoid importing all infos into a list, which is ram-wasting
		import linecache
			print(linecache.getline(r'D:\z.txt',10))
			# line cache could result in large cache space
			
		def getline(thefilepath, desired_line_number):
			if desired_line_number < 1: return ''
			for current_line_number, line in enumerate(open(thefilepath, 'rU')):
				if current_line_number == desired_line_number - 1 : return line
			return ''
			
			
		while ind.isalpha():
			# set pannel to be 500bp
			# as long as smaller than max ion read length
			ind=fafh.read(500)
			numind += 1
			fafh.seek(numind,0)
			dctfh.write(ind+"\n")
			
			#note that for seek(offset,wherence) wherence=0 from that absolute beginning, 1 from current place, 2 from ending
	
0722
	###################################
	# re compilation begin here #######
	softre=re.compile(r'.*?(\d+)S.*?')											# soft number
	insre=re.compile(r'.*?(\d+)I.*?')											# insertion number
	matre=re.compile(r'.*?(\d+)M.*?')											# match number
	padre=re.compile(r'.*?(\d+)P.*?')											# padding number
	mismre=re.compile(r'.*?(\d+)X.*?')											# mismatch number
	mat2re=re.compile(r'.*?(\d+)\=.*?')											# type 2 match number	
	###################################
	
	tmap index -vf ucsc.hg19.chrX.std.fasta
	tmap map1 -f /path/to/ref.fasta -r candidate.fa -s output.sam
	tmap mapall -f ucsc.hg19.fasta -r test.fa -s test.mapall.sam stage1 map4
	
	
	
0724
	# findrobots.py
	import gzip
	import io
	import glob
	def find_robots(filename):
		'''
		Find all of the hosts that access robots.txt in a single log file
		'''
		robots = set()
		with gzip.open(filename) as f:
			for line in io.TextIOWrapper(f,encoding='ascii'):
				fields = line.split()
				if fields[6] == '/robots.txt':
					robots.add(fields[0])
		return robots
		
	# normal version
	def find_all_robots(logdir):
		'''
		Find all hosts across and entire sequence of files
		'''
		files = glob.glob(logdir+'/*.log.gz')
		all_robots = set()
		for robots in map(find_robots, files):
			all_robots.update(robots)
		return all_robots
	
	if __name__ == '__main__':
		robots = find_all_robots('logs')
		for ipaddr in robots:
			print(ipaddr)
			
	# comcurrent.futures
	from concurrent import futures
	def find_all_robots(logdir):
		'''
		Find all hosts across and entire sequence of files
		'''
		files = glob.glob(logdir+'/*.log.gz')
		all_robots = set()
		with futures.ProcessPoolExecutor() as pool:
			for robots in pool.map(find_robots, files):
				all_robots.update(robots)
		return all_robots
	
	# Processing pool (see below for initiazation)
	pool = None
	# Performs a large calculation (CPU bound)
	def some_work(args):
		...
		return result
	# A thread that calls the above function
	def some_thread():
		while True:
			...
			r = pool.apply(some_work, (args))
			...
	# Initiaze the pool
	if __name__ == '__main__':
		import multiprocessing
		pool = multiprocessing.Pool()
	
	
0725
	testing concurrent.futures.ProcessPoolExecutor
	obtain results from mapped.sorted.bams
		
	test@218.19.14.181
	ig34438011
	
0726
	from multiprocessing import Pool as ThreadPool								# multiprocessing core
	from multiprocessing.dummy import Pool as ThreadPool						# multiple threads
	from concurrent.futures import ProcessPoolExecutor as ThreadPool			# multiple threads
	
	for python strip, only changes the output, while leaving ori line impact
	
	for tmap, use -a to specify output mode
	
	
	DMD:
		 NM_004006       chrX    -       31137344        33229673        31140035        33229429        79
	
	
	ps -elf | grep 'ljzhang' | grep 'python' | awk '{print $4}' | xargs kill -9
	
	do not use range() or in list operation, for this is time consuming, use compare (><=) instead
	
	print only takes string, while file reading takes strings
	finished filtra.py and get_breakpoint.py
	
	ipython installed
	ipython
	
	
0728
	get samtools reads enrichment infos
	get histplot out for different xlab axis
	> paste(dep3[1,1],dep3[1,2],sep='-')
	[1] "31137232-31138232"
	
	strsplit(<str>,split=<sep>)
	paste(<str>…, sep = “ “)
	
	barplot(beside=TRUE,dep3$depth,dep3$intv,border='black',space=rep(2,nrow(dep3)),width=rep(2,nrow(dep3)),cex.axis=.5,ylim=c(0,max(dep3[,3])))

	bp <- ggplot(你的图像之类)     
	bp + theme(axis.text.x = element_text(angle = 70, hjust = 0.5, vjust = 0.5))
	
	/home/ljzhang/project/dmd/0707_bam_test/test_gatk_pipe_TMP/breakpointTest
	retmapped.reshaped->filtra.py->.anno->get_breakpoint_info.py->.anno.breakinfo info with soft-clipping position included
	/home/ljzhang/project/dmd/0707_bam_test/test_gatk_pipe_TMP/dmdDepth:
	tmapped.sorted.bam->samtools_depth.py->dmddepth->depth_interval.py->gap.txt->dmd_breakpoint_plot.R
	
	
0803
	using R file and Sys module for file and dir operation
	getting plot out for +- gap*times
	testing with 100kb gap threshold for breakpoint
		results are good
		now start merge and print out line plots
		
	use:
		a<<-3 or assign("a",3,envir=.GlobalEnv)
	reconstruct procedure of a dmd paper by BIG 
	install cmd: requiring sudo
		sudo apt-get install libcairo2-dev
		sudo apt-get install libxt-dev
		
	note that chromosome positions are known to be not continuous after smatools depth
	
	note that in 1k window the diff function would not work as a classifier
	another problem is that the breakpoint maybe startpoint sensitive
	
	try to optimize the 1k points version instead of the 1k gap
	
	
0805
	after normalization is performed, a z score estimation is applied for breakpoint search
	
	
0808
	first finish using zscore for breakpoint-sniffer
	then merge with python pipeline
	results are needed!!!!!!!!!!
	be responsible for your script
	
	R ave function
	
0809
	kNN classifier number not defined and hard to reach with only present position and 1k position downstream
	one posible measure is to add CIGAR string 
	
	ISODATA algorithsm obtains a relative resonable k by merging and splitting data clusters
	
	computing time consumption
	
	RSICNV: https://github.com/yhwu/rsicnv
	installed in 45 server: auto-GC-adjustment
	
	a modified work flow:
	1. generate 1*gap positions on chrX DMD target gene area, include a times*gap position info on tail
	2. use biend search method, search for breakpoint using tail dataframe info
	3. apply a k-means classification
		1> use zscore for a dummy coding to facilate filtration
		2> dummy coding requires infos including cv mean and zscore, use univ-info as control
	
	
0810
	kmeans on line, > 9h
	reconsider 1k window, use dummy coding for missing values (complete all positions)
	prof:
		1. only 238 postions sequenced between subset(i3,b>31780000 & b<31800000)
		2. first calculate 2-point interval for overall data
		3. understand procedure for GC-adjustment
	only problem is that 1k amy contain lots of indels, thus use interval between 2 start point to filtration
	
	
	for ionadmin data auto-processing:
		variation priority: deletion & point mutation > duplications
		merge continuous info
		homzygosity info are stored in both annotation file and vcf file in tvcinfo folder
		script should be rewritten in info_sum
	
	the complete 1k window procedure is most time consumming
	
	gc adjustment is used for CNV detection
	
	one possible solution to reduce kmeans time is to apply classification for each intervals
	
	calculate the between point length
		1> under hypothesis that at deletion place, reads coverage are more sparse
		2> calculate and reads and sparse divergence
		3> use python
		
0811
	use blast for reads alignment
	
	use mapping quality as filtration control
	
	use CIGAR string for score measuring
	numba and numbapro is used for CUDA parallelization
	
	
	
0815
	Cairo package works for a single window pannel, thus put library(Cairo) inside of a target function
	
0816
	blast packages installed
	
0817
	get prepared for the monthly report, 
	
	
0819
	get prepared for the paper review
	automize pipeline from raw data to semiautomized results
	
	
0823
	THetA: inferring intra-tumor heterogeneity from high-throughput DNA sequencing data:
		http://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-7-r80
	
	data stored on 47 server: /home/ljzhang/project/dmd/0725
		cmd: scp *rawlib.basecaller.bam ljzhang@192.168.0.45:/home/ljzhang/project/dmd/0707_bam_test
	
0825
	solving SUP documentations of breakpoint_caller
	test dmd_breakpoint_caller_0819.py
	caller program is on run now
	
	finish that SUP docu
	
0826
	finish a SUP ocument 
	prepare for the paper presentation
	
0829
	finished a dmd_breakpoint_caller_0819.py SOP document
	new dmd batch sequencing results are available at:
		35 ionadmin:
			/results/analysis/output/Home/Auto_sn247560054_SN2-409-P35-HiQ-DMD-10za1-Pooling-160828_608_942/basecaller_results
			'barcode6是22exon重复，barcode11是3-17重复，barcode8、9、10是正常男性对照，barcode12、13、14是正常女性，57是陈婷全血，58是陈婷羊水'

	a automatic reporter script is put at:
		45 server: /home/ljzhang/project/35_server_dmd_report.py
	note that the great variation of dmd depth can not be corrected by normal data, for this is experiment biased
	
	method cange: 
		1 change our plotting method to refine a 20% percent range-span 
		2 change out line plotting to histograms
	
	
0830
	use python to rewrite the mean searching method for deletion breakpoint
	change ploting method to barplot, use examples in TCGA patient pattern
	
	SOFTWARES FOR DETECTION CNV:
	Translocation and Inversions:	
		Discordant paired end
			BreakDancer
			Hydra
			VariationHunter
			PEMer
			GASVPro
			
			CREST
			Slope
	
	CNV
		Raw coverage based
			EWT
		Coverage rato based
			SeqSeq
			CNVnator
			CNAseg
			CNV-seq
			CONTRA
			CoNVEX
			ExomeCNV
	
	Insertions and Deletions
			Pindel
			GATK
			VarScan2
			Dindel
			Stampy
			SAMtools
		
	CNV
		Raw coverage based
	EWT Whole genome only, does not require normal control http://rdxplorer.sourceforge.net
		Coverage ratio based
	SeqSeq Whole genome only; requires normal control http://www.broadinstitute.org/software/cprg/?q=node/39
	CNVnator Whole genome only; requires normal control http://sv.gersteinlab.org
	CNAseg Whole genome only; requires normal control http://www.compbio.group.cam.ac.uk/software/cnaseg
	CNV-seq Whole genome only; requires normal control http://tiger.dbs.nus.edu.sg/cnv-seq/
	CONTRA Exome or targetd panels; requires normal controls http://sourceforge.net/projects/contra-cnv/
	CoNVEX Exome; requires normal controls ftp://ftp.sanger.ac.uk/pub/users/pv1/CoNVex/Docs/CoNVex.pdf
	ExomeCNV
		Exome or targeted panels; requires normal controls; evaluates B-allele frequency 
		https://secure.genome.ucla.edu/index.php/ExomeCNV_User_Guide
	
	EWT: event wise testing
		RDXplorer
		
		should be mapped sorted and samtools rmdup
			samtools rmdup IonXpress_006_rawlib.basecaller.sorted.bam out.rdup.bam
		
0831:
	added soft-clipping beginning to breakinfo_filtra_0831_test.py
	excutable
	
	to avoid the samtools time cosuming step
	use dmddepth file and conducting a bi-fork searching and store a dict
	
	
	use python traceback module for tracing errors, better than except Exception as error: print error
		
		import sys
		import traceback
		import test1
		
		a=10
		b=0
		
		try:
			print test1.division(a,b)
		except:
			print 'invoking division failed.'
			traceback.print_exc()
			sys.exit(1)
			
0901:
	note that python passes its params by soft-reference
	just calm down for writing python functions
	
	modify potential point filtration for V1.02, for considerring head soft-clipping
	
	use sorting and classification method for breakpoint filtration
	consider direction of mapping using sam file info, use this as a filtration
	
	use 500bp as a classification threshold
	
	make notes of samples processed, extract info from DMD PLUG-IN, keep records for future analysis
	
0905:
	reconsider soft-clipping and FLAG info
	if softclipping length<20 or raw-reads<20:
		filtration
	
	if rawreads mapping FLAG != raw-reads cut mapping FLAG:
		filtration
	
	put annotation afterwards
	
	cutting python file afterwards
	
0906:
	for FLAG 16:
		still treat giving reads as positive string
	
	
	
0907:
	breakpoint_filtra_0901_v102.R has barplot method
	avoid hist with big array
	
	use color and filtration
	
	
0908
	col added using an extra col 5,8 indicates pink and grey
	there would be monthly report
	careful with putting if after an output filehandle
	use traceback module for tracking error info
	# py
	import traceback
	traceback.print_exc()
		
		http://blog.csdn.net/csujiangyu/article/details/45332083
	
	a useful tool for filtra out replicates for diverse data_struct types
	
	def dedupe(items,key=None):
		seen=set()
		for item in items:
			val=item if key is None else key(item)
			if val not in seen:
				yield item
				seen.add(val)
				
	
	from operator import itemgetter
	rows_by_fname=sorted(rows,key=itemgetter('fname'))
	rows_by_fname=sorted(rows,key=lambda r:r['fname'])
	## note that itemgetter could accept multiple keys for sorting
	rows_by_lfname=sorted(rows,key=itemgetter('lname','fname'))
	
	from operator import attrgetter
	which has similar function as itemgetter
	
	from operator import itemgetter
	from itertools import groupby
	
	rows.sort(key=itemgetter('data'))
	
	for date,items in groupby(row,key=itemgetter('data'))
		print(date)
		for i in items:
			print(' ',i)
	
	
0909
	v1.05 seems not working well with 07 deletion samples
	
	careful with using classification with deletion breakpoint detection
	
	
	
	awk '$4>31725991 && $4<31779286{print $0}' ../IonXpress_002_rawlib.basecaller.chrX.sorted.sam | wc -l
	200
	awk '$4>31725991 && $4<31779286{print $0}' IonXpress_002_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | wc -l
	25
	
	awk '$4>31754549 && $4<31954212{print $0}' ../IonXpress_003_rawlib.basecaller.chrX.sorted.sam | wc -l
	226
	awk '$4>31754549 && $4<31954212{print $0}' IonXpress_003_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | wc -l
	47
	
	awk '$4>31810925 && $4<31964877{print $0}' ../IonXpress_004_rawlib.basecaller.chrX.sorted.sam | wc -l
	743
	awk '$4>31810925 && $4<31964877{print $0}' IonXpress_004_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | wc -l
	198
	
	awk '$4>31803171 && $4<31907486{print $0}' ../IonXpress_005_rawlib.basecaller.chrX.sorted.sam | wc -l
	815
	awk '$4>31803171 && $4<31907486{print $0}' IonXpress_005_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | wc -l
	150
	
	awk '$4>32398857 && $4<32914837{print $0}' ../IonXpress_007_rawlib.basecaller.chrX.sorted.sam | wc -l
	745
	awk '$4>32398857 && $4<32914837{print $0}' IonXpress_007_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | wc -l
	47
	
	awk '$4>31698235 && $4<31926191{print $0}' ../IonXpress_008_rawlib.basecaller.chrX.sorted.sam | wc -l
	343
	awk '$4>31698235 && $4<31926191{print $0}' IonXpress_008_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | wc -l
	93
	
	great decrease of supporting reads when detecting large deletion
	
	
	
	
0914:
	perl -le '$a="abc";$a=~tr/a-z/A-Z/;print $a'
	try to cut chrX dmd gene out to make a reference with 5kb up-down stream
	
	samtools faidx ucsc.hg19.fasta chrX:31132344-33234673 > ucsc.hg19.dmd.5kb.fasta
	
	for 35 server:
		python 35_server_dmd_report.py -D /results/analysis/output/Home/sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_re_982/plugin_out/DMD_plugin_out.5852 -O 35_lj_0914_dmd_results.txt
		should wait until processing is over
	
	awk '$4>31810625 && $4<31811625{print $4,$6,$21,$23,$26,$28}' OFS='\t' IonXpress_004_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped | less -S
	
	
0920
	prepare for rewriting DMD-plugin
	use negative sample as control
	prepare model using Normalization filtration
	
	
	
	
0921
	sum up method abstract for CNV search algorithsm
	make a ppt
	rethink about DMD PLUGIN, make pre for new PLUGIN searching method
	
	Pindel:
		git clone git://github.com/genome/pindel.git
		cd pindel 
		./INSTALL /path/to/samtools_FOLDER/
		# test the code
		cd demo
		../pindel -f simulated_reference.fa -i simulated_config.txt -o output
	
	cnD(optional):
		sudo apt-get install gcd
		sudo yum install tango-devel
		sudo cpan Set::IntSpan
		
		cd /programinstallers/
		wget -N ftp://ftp.sanger.ac.uk/pub4/resources/software/cnd/cnD-1.3.tar.gz
		tar -zxvf cnD-1.3.tar.gz
		
	
	
	pindel works for ./demo examples
	help page here:
	http://gmt.genome.wustl.edu/packages/pindel/quick-start.html
	
	pindel seems to only works with small deletion?
	
	error happend when applying bam2pindel.pl for getting pindel format txt
		bam2pindel.pl -i IonXpress_007_rawlib.basecaller.sorted.bam -o ./out/prefix -s ion007 -pi 250 -om
		no reads writed for torrent data
	
	
0922
	using python for getting pipeline output
	
		import subprocess as sp
		import sys

		pipe=sp.Popen('samtools view IonXpress_001_rawlib.basecaller.sorted.bam', stdout=sp.PIPE,shell=True)
		ind=0
		outfh=open('test.txt','w')
		for line in iter(pipe.stdout.readline,''):
			print(line.strip())
			outfh.write(line)
			ind+=1
			if ind>20:
					sys.exit(1)
	
	
	use Q30 as a threshold for reads filtration
	
	GC-value was calculated according to reference genome
	
	use index and bifork searching
	
	
	to get a GC content from reference genome, file.seek method is required
	
		def generateDict(chrfile,outfile):
			fafh=open(chrxfile,'r')
			dctfh=open(dictfile,'w')
			
			# generating dict #
			ind=''
			ind=fafh.read(1)
			fafh.seek(0,0)
			
			numind=0
			
			# putting all pannel in one file, if reading in is necessary, use readlines()
			while ind.isalpha():
				# set pannel to be 500bp
				# as long as smaller than max ion read length
				ind=fafh.read(500)
				numind += 1
				fafh.seek(numind,0)
				dctfh.write(ind+"\n")
			fafh.close()
			dctfh.close()
	
	
0923:
	prepare well for DMD project
	
	DMD plug-in re-run
		Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615
	
	plug-in in queue
	original file found at:
		D:\thor\documents\45capitalgenomics\dmdProject\DMD项目一二批样本\陈雅莉20160831-DMD&KD测序样品信息_0923_
	
	cancer heterogeneity and aggressive cancer
	
	
0926:
	keep track with DMD data
	
	list().index(element)
	[].extend([])
	
	finished a EWT algorhism
	
	
0927:
	problem with EWT comes from sample
	first of all: GC adjustment may not perform well with 0 depth, maybe because of targeted area
	maybe we should exclude no depth region out from calculating median depth
	
	ll -hl Ion*.bam | awk '$5~/G/ {print $8}'
	
	move 20160831_615_955 sample to 45 for analysis
	
	questions about applying EWT algorithsm:
		keep all 100bp window under a normal distribution
		whether DMD region is targeted for other region depth is near 0
		DMD region span chrX:31132344-33234673
		answer is yes, only DMD region and a small region of chrX is designed
		
	DMD whole gene bed noted on 35server: X       31137345        33229673
	
	we rerun 0831_615_955 data with breakpoint caller, figures show that this dataset is based on whole exon sequencing which is not suitable for breakpoint sniffer
	
	
	
0928:
	use simulating depth info to test the power of EWT algorithsm
	note that a specific region of chrX is selected as an inner reference for detecting CNV in DMD
	
	use loc in pandas to avoid the copying of a slice
	df.loc[df['pos']>10,('pos')]=0
	df.loc[(df['pos']>10) & (df['pos']<20),('pos')]=0
	
	pandas is awesome
	
	
0929:
	format=lambda x: '%.2f' % x
	frame.applymap(format)
	
	frame.sort_index()
	
	frame.sort_index(by=['a','b'])
	# error produced in Ipython, recommanded for using 
	frame.sort_values(by='zsc',ascending=False)
	
	a key factor is to mark PCR duplications on 
	
	$. in perl returns current row nbr, started from 1
	
	use picard markduplicate for remove duplicates
	samtools rmdup also could be used for remove duplicates
	
	############################
	## python pandas packages ##
	############################
	
	for Series:
	
		obj=Series([])
		obj.values
		obj.index
		
		Series([],index=[])
		
		'b' in obj
		
		sdata={}
		obj3=Series(sdata)
		
		pd.isnull(obj)
		pd.notnull(obj)
		obj.isnull()
		
		obj.name=''
		obj.index.name=''
		
	for DataFrame:
		
		DataFrame(data,columns=[],index=[])
		
		frame['abc']
		frame.abc
		frame.ix['']
		
		frame['extern']=frame.state=='Ohio'
		
		frame.T
		
		obj.reindex([],fill_value=0)
		
		frame=DataFrame(np.arange(9).reshape((3,3)),index=['a','c','d'],cilumns=['Ohio','Texas','California'])
		
		frame.reindex(columns=states)
		
		frame.ix[[],states]
		
		obj.drop(['d','e'])
		obj['b':'c']=5
		
		data.ix[[],[]]
		df1.add(df2,fill_value=0)
		
		f=lambda x:x.max()-x.min()
		frame.apply(f)
		frame.apply(f,axis=1)
	
		def f(x):
			return Series([x.min(),x.max()],index=['min','max'])
			
		frame.apply(f)
		
		format=lambda x:'%.2f' % x
		frame.applymap(format)
		
		frame['e'].map(format)
		
		obj=Series(range(4),index=[])
		obj.sort_index()
		
		frame.sort_index(axis=1)
		frame.sort_index(by=['a','b'])
		
		obj.rank(method='first')
		obj.rank(ascending=False,method='max')
		
		df.sum()
		df.sum(axis=1)			# for rows
		df.mean(axis=1,skipna=False)
		
		df.idxmax()
		df.cumsum()
		
			subfunctions:
			count
			describe
			min,max
			argmin,argmax
			idxmin,idxmax
			quantile
			sum
			mean
			median
			mad
			var
			std
			skew
			kurt
			cumsum
			cummin,cummax
			cumprod
			diff
			pct_change
		
		returns.MSFT.corr(returns.IBM)
		
		returns.MSFT.cov(returns.IBM)
		
		returns.corr()
		returns.cov()
		
		returns.corrwith(returns.IBM)
		
		uniques=obj.unique()
		obj.value_counts()
		
		pd.value_counts(obj.values,sort=False)
		mask=obj.isin(['b','c'])
		
			dropna
			fillna
			isnull
			notnull
			
		pd.merge(df1,df2,on='key',how='left')
		
		
0930:
	length of all hg19 chromosomes
		chr1	249250621
		chr2	243199373
		chr3	198022430
		chr4	191154276
		chr5	180915260
		chr6	171115067
		chr7	159138663
		chrX	155270560
		chr8	146364022
		chr9	141213431
		chr10	135534747
		chr11	135006516
		chr12	133851895
		chr13	115169878
		chr14	107349540
		chr15	102531392
		chr16	90354753
		chr17	81195210
		chr18	78077248
		chr20	63025520
		chrY	59373566
		chr19	59128983
		chr22	51304566
		chr21	48129895
		chr6_ssto_hap7	4928567
		chr6_mcf_hap5	4833398
		chr6_cox_hap2	4795371
		chr6_mann_hap4	4683263
		chr6_apd_hap1	4622290
		chr6_qbl_hap6	4611984
		chr6_dbb_hap3	4610396
		chr17_ctg5_hap1	1680828
		chr4_ctg9_hap1	590426
		chr1_gl000192_random	547496
		chrUn_gl000225	211173
		chr4_gl000194_random	191469
		chr4_gl000193_random	189789
		chr9_gl000200_random	187035
		chrUn_gl000222	186861
		chrUn_gl000212	186858
		chr7_gl000195_random	182896
		chrUn_gl000223	180455
		chrUn_gl000224	179693
		chrUn_gl000219	179198
		chr17_gl000205_random	174588
		chrUn_gl000215	172545
		chrUn_gl000216	172294
		chrUn_gl000217	172149
		chr9_gl000199_random	169874
		chrUn_gl000211	166566
		chrUn_gl000213	164239
		chrUn_gl000220	161802
		chrUn_gl000218	161147
		chr19_gl000209_random	159169
		chrUn_gl000221	155397
		chrUn_gl000214	137718
		chrUn_gl000228	129120
		chrUn_gl000227	128374
		chr1_gl000191_random	106433
		chr19_gl000208_random	92689
		chr9_gl000198_random	90085
		chr17_gl000204_random	81310
		chrUn_gl000233	45941
		chrUn_gl000237	45867
		chrUn_gl000230	43691
		chrUn_gl000242	43523
		chrUn_gl000243	43341
		chrUn_gl000241	42152
		chrUn_gl000236	41934
		chrUn_gl000240	41933
		chr17_gl000206_random	41001
		chrUn_gl000232	40652
		chrUn_gl000234	40531
		chr11_gl000202_random	40103
		chrUn_gl000238	39939
		chrUn_gl000244	39929
		chrUn_gl000248	39786
		chr8_gl000196_random	38914
		chrUn_gl000249	38502
		chrUn_gl000246	38154
		chr17_gl000203_random	37498
		chr8_gl000197_random	37175
		chrUn_gl000245	36651
		chrUn_gl000247	36422
		chr9_gl000201_random	36148
		chrUn_gl000235	34474
		chrUn_gl000239	33824
		chr21_gl000210_random	27682
		chrUn_gl000231	27386
		chrUn_gl000229	19913
		chrM	16571
		chrUn_gl000226	15008
		chr18_gl000207_random	4262
	
	
	
	
1008:
	df[-1:].Pos
	stats.norm.ppf(1-0.000001)
	np.concatenate([arr,arr],axis=1)
	pd.concat([s1,s2,s3],axis=0)
	
	samtools depth IonXpress_003_rawlib.basecaller.sorted.bam -r chrX:33229684 | awk '$3>80' | less -S
	
	a big problem is that when applying inner reference data, sequencing range is not suffisant
	also the ref region for different sample differs
	
	samtools depth IonXpress_009_rawlib.basecaller.sorted.bam -r chrX:33229684 | awk '$3>40' | less -S
	
	
1009:
	get potential ref Region for each sample
	
	DMD results are ready for analysis:
		37 server:
			8 sample for DMD breakpoint
			/results/analysis/output/Home/Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018
			16 clinical sample
			/results/analysis/output/Home/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_tn_1011
			
	
	
1010:
	import time
	time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))
		return str result
	
1011:
	json file format seems to be wrong
	modify perl reading or rewrite one yourself
	
	rewrite pick_out_dmd_bams_from_a_run.py
	getinto pipeline
	
	
1012:
	python dmd_plugin_pipeline.py -P /media/disk2/ljzhang/project/DMD_plugin/DMD_plugin_test -R /media/disk2/ljzhang/project/DMD_plugin/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/plugin_out/DMD_plugin_out.1530 -B /media/disk2/ljzhang/project/DMD_plugin/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/basecaller_results -A /media/disk2/ljzhang/project/DMD_plugin/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010
	
	
	
1013:
	Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]
	http://cancer.sanger.ac.uk/cosmic/search?q=DLBCL#samp
		
	python dmd_plugin_pipeline.py -P /media/disk2/ljzhang/project/DMD_plugin/DMD_plugin_test -R /media/disk2/ljzhang/project/DMD_plugin/Auto_user_1013/plugin_out/DMD_plugin_out.1530 -B /media/disk2/ljzhang/project/DMD_plugin/Auto_user_1013/basecaller_results -A /media/disk2/ljzhang/project/DMD_plugin/Auto_user_1013
	
	
	
	cd /media/disk2/ljzhang/project/DMD_plugin/Auto_user_1013
	
	../DMD_plugin_test/bam2fastq-1.1.0/bam2fastq ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib.bam ./plugin_out/DMD_plugin_out.1530/IonXpress_006_rawlib.fastq -q
	
	../DMD_plugin_test/FastQC/fastqc -f fastq -o ./plugin_out/DMD_plugin_out.1530/fastqc_result -t 1 ./plugin_out/DMD_plugin_out.1530/IonXpress_006_rawlib.fastq
	
	rm ./plugin_out/DMD_plugin_out.1530/IonXpress_006_rawlib.fastq
	
	# perl ${PLUGIN_PATH}/BamStat.v3.pl -regionFile ${PLUGIN_PATH}/DMD100_Targets.bed -bamFile $bam -outDir ${RESULTS_DIR}/bam_stat_first -plot -largeDeletion -largeInsertion
	perl ../DMD_plugin_test/BamStat.v3.pl -regionFile ../DMD_plugin_test/DMD100_Targets.bed -bamFile ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib.bam -outDir ./plugin_out/DMD_plugin_out.1530/bam_stat_first -plot -largeDeletion -largeInsertion
	
	# perl ${PLUGIN_PATH}/annovar/table_annovar.pl -protocol refGene -buildver hg19 -operation g -nastring - -remove ${RESULTS_DIR}/bam_stat_first/${prefix}_large_deletion.txt -outfile ${RESULTS_DIR}/bam_stat_first/${prefix}_large_deletion_annotation.txt ${PLUGIN_PATH}/annovar/humandb/
	perl ../DMD_plugin_test/annovar/table_annovar.pl -protocol refGene -buildver hg19 -operation g -nastring - -remove ./plugin_out/DMD_plugin_out.1530/bam_stat_first/IonXpress_006_rawlib_large_deletion.txt -outfile ./plugin_out/DMD_plugin_out.1530/bam_stat_first/IonXpress_006_rawlib_large_deletion_annotation.txt ../DMD_plugin_test/annovar/humandb/
	# no success
	
	# perl ${PLUGIN_PATH}/BamStat.v3.a.2015.5.29.pl -regionFile ${PLUGIN_PATH}/DMD.enst7033.exon.bed -bamFile $bam -outDir ${RESULTS_DIR}/bam_stat_second -plot -gcCorrection
	perl ../DMD_plugin_test/BamStat.v3.a.2015.5.29.pl -regionFile ./DMD_plugin_test/DMD.enst7033.exon.bed -bamFile ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib.bam -outDir ./plugin_out/DMD_plugin_out.1530/bam_stat_second -plot -gcCorrection
	# Smartmatch is experimental at ../DMD_plugin_test/BamStat.v3.a.2015.5.29.pl line 454.
	# No such file or directory at ../DMD_plugin_test/BamStat.v3.a.2015.5.29.pl line 72.

	
	# perl ${PLUGIN_PATH}/BamStat.v3.flank.pl -regionFile ${PLUGIN_PATH}/DMD.enst7033.exon.including.flanking.10bp.bed -bamFile $bam -outDir ${RESULTS_DIR}/bam_stat_third -plot -flank
	perl ../DMD_plugin_test/BamStat.v3.flank.pl -regionFile ../DMD_plugin_test/DMD.enst7033.exon.including.flanking.10bp.bed -bamFile ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib.bam -outDir ./plugin_out/DMD_plugin_out.1530/bam_stat_third -plot -flank
	# file in third, but reports several perl syntex errors
	
	
	# perl ${PLUGIN_PATH}/BamStat.v3.a.pl -regionFile ${PLUGIN_PATH}/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed -regionFileDmd ${PLUGIN_PATH}/dmd_region.bed -regionFileRef ${PLUGIN_PATH}/ref_region.bed -bamFile $bam -outDir ${RESULTS_DIR}/bam_stat_fourth_a -gcCorrection -sex m
	perl ../DMD_plugin_test/BamStat.v3.a.pl -regionFile ../DMD_plugin_test/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed -regionFileDmd ../DMD_plugin_test/dmd_region.bed -regionFileRef ../DMD_plugin_test/ref_region.bed -bamFile ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib.bam -outDir ./plugin_out/DMD_plugin_out.1530/bam_stat_fourth_a -gcCorrection -sex m
	
	
	### CNV calling here
	# perl ${PLUGIN_PATH}/calculate_zscore_intra_run_v2.pl ${RESULTS_DIR}/bam_stat_fourth_a ${RESULTS_DIR}/dmd_project/barcode2sample.txt ${PLUGIN_PATH}/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed ${RESULTS_DIR}/CNV_result
	perl ../DMD_plugin_test/calculate_zscore_intra_run_v2.pl ./plugin_out/DMD_plugin_out.1530/bam_stat_fourth_a ./plugin_out/DMD_plugin_out.1530/dmd_project/barcode2sample.txt ../DMD_plugin_test/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed ./plugin_out/DMD_plugin_out.1530/CNV_result
	# throws error: Can't stat /usr/local/lib/site_perl
	
	
	############# barcoding mkdirs #############
	
	
	
	# mutation annotation starts here
	# $PLUGIN_PATH/samtools sort -p 2 $bam ${RESULTS_DIR}/dmd_project/${bam_prefix}_sorted
	# $PLUGIN_PATH/samtools index ${RESULTS_DIR}/dmd_project/${bam_prefix}_sorted.bam
	# $PLUGIN_PATH/tvc --parameters-file ${PLUGIN_PATH}/variant_caller_scripts/json/startplugin.mlld2.from.second.tvc4.for.hom.json --reference /results/referenceLibrary/tmap-f3/hg19/hg19.fasta --input-bam ${RESULTS_DIR}/dmd_project/${bam_prefix}_sorted.bam --target-file ${PLUGIN_PATH}/DMD100_Targets.extend.200bp.flanking.chr.bed --output-dir ${RESULTS_DIR}/tvc_out --output-vcf $ionXpress.vcf
	
	../DMD_plugin_test/samtools sort -p 2 ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib.bam ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib_sorted
	
	../DMD_plugin_test/samtools index ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib_sorted.bam
	
	
	######################### important tips #########################
	### when using tvc, export whole tvc folder
	# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PLUGIN_PATH/tvc_lib
	
	../DMD_plugin_test/tvc --parameters-file ../DMD_plugin_test/variant_caller_scripts/json/startplugin.mlld2.from.second.tvc4.for.hom.json --reference /results/referenceLibrary/tmap-f3/hg19/hg19.fasta --input-bam ./plugin_out/DMD_plugin_out.1530/dmd_project/IonXpress_006_rawlib_sorted.bam --target-file ../DMD_plugin_test/DMD100_Targets.extend.200bp.flanking.chr.bed --output-dir ./plugin_out/DMD_plugin_out.1530/tvc_out --output-vcf IonXpress_006.vcf
	
	
	
	
1014:
	test all commands on 37 server
	
	shame that the error was caused because of using perl instead of python
	
	
	
1017:
	bed file: DMD85_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed
	
1018:
	prepare for rewriting calculating z-score part
	
	several things that should be paid attention to:
		keep in close contact with YL, keeping tracking with DMD clinical process
		keep everything simple and slim
		dont ever be that moron again, get yourself a plan for your work
		
	for intra-run zscore calculating script:
		1\ Reading in bed file
		2\ Mapping and process raw genome with BED positions
		3\ get GC spectrum for whole DMD region
		4\ linearize GC content-depth relationship
		5\ use q20 as threshold for summerize depth info
		6\ calculating realtive depth to ref region
		7\ intra run relative depth, mean and std
		8\ calculating z-score
		9\ redo using GC corrected data
		
	correct string to upper case for calculating GC content in EWT SCRIPT
	
	
1019:
	prepare for DMD history file
		Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630										no results(in run with older version			finished at 2016 1020)
			35:/results/analysis/output/Home/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980
		2016_09_12
		
		Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_tn					no results(in run with older version			finished at 2016 1020)
			35:/results/analysis/output/Home/Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955
		2016_08_31
			
		Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_tn						results ready
			37:/results/analysis/output/Home/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010
		2016_09_30
		
		Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621						results ready
			37:/results/analysis/output/Home/Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018
		2016_10_08
		
		Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_tn								results ready
			35:/results/analysis/output/Home/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037
		2016_10_14
		
		Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_tn_377/							results ready
			34:/results/analysis/output/Home/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376
		2016_08_25
		
1020:
	vt[(vt.index>2) & (vt.index<5)]
	
	summary for pd:
		Series:
			.index
			pd.isnull()
			pd.notnull()
			.reindex([],method='ffill/pd OR bfill/backfill')
				fill_value
				limit
				level
				copy
			.drop([])
			.sort_index()
			.rank()
			.rank(method='first')
			.count
			.describe
			.min
			.max
			.argmin
			.argmax
			.idxmin
			.idxmax
			.quantile
			.sum
			.mean
			.median
			.mad
			.var
			.std
			.skew
			.kurt
			.cumsum
			.cummin
			.cummax
			.cumprod
			.diff
			.pct_change
			.unique
			.value_counts
			.isin([])
			.dropna()
				
		index:
			append
			diff
			intersection
			union
			isindelete
			dropinsert
			is_monotonic
			is_unique
			unique
			
		DataFrame:
			.ix[]
			['']
			.T
			.add(df2,fill_value=0)
			.sub/div/mul
			DataFrame(,columns=[],index=[])
			.apply(,axis=)	# 	1 for row, default to be applied to columns
			.applymap(format)
			.sort_index(axis=1,ascending=False)
			.sort_index(by='b')
				by=[]
			
		pd.merge(df1,df2,on='key')
		pd.merge(df1,df2,left_on='lkey',right_on='rkey')
		pd.merge(df1,df2,how='outer/inner/left')
		
			
	note that multiple object could be returned from a python function, in form of a tuple
	
	use relative depth for filtration of large indel of duplications
	
	
1024:
	note that np.round([],digits)
	
	
1026:
	Series.sort_index()
	Series.sort_values(ascending=False)
	
	pd.DataFrame.sort_index()
	pd.DataFrame.sort_index(by='')
	
1031:
	get cnv preparing data out and prepare report for presenting
	searching for way for IDE connecting servers
	one thing about upgrades version is that:
		should be sure that the older version has no syntax or logical error
		
	in future versions, correct depth info using samtools depth
	
	in a recent DMD run, reference is selected, this change has beem corrected
	
	selecting samples fr 2M pannel
	
	for multiple annotation case:
		avoid unpredictable problem by writing full path
		
		
1101:
	get test data ready for inspecting
	get total trends out
	
	Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697
	
	/results/analysis/output/Home/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/plugin_out/DMD_plugin_out.6040/tvc_out/IonXpress_080
	
	/results/analysis/output/Home/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/plugin_out/DMD_plugin_out.6040/tvc_out/IonXpress_079
	
	/results/analysis/output/Home/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/plugin_out/DMD_plugin_out.6009/tvc_out/IonXpress_012
		
		
	for looped functions in python:
		use all params passing to avoid global syntex error
		
	x=10
	def f1():
		x=7
		print x
		f2()
		
	def f2():
		print x+3
		
	f1()

	this returns 7,13
	
	x=10
	def f1():
		x=7
		print x
		f2(x)
		
	def f2(x1):
		print x1+3
		
	f1()
	
	this returns 7,10
	
	#############################################
	# getting depth info and relDepth info plot 
	#############################################
	
	setwd('D:\\thor\\documents\\45capitalgenomics\\dmdProject\\new_DMD_plugin\\data_test')
	
	refDmd=read.table('refDmdInfo.txt',header=T,stringsAsFactors=F)
	univDmd=read.table('univDmdInfo.txt',header=T,stringsAsFactors=F,fill=NA)
	unvGc=read.table('unvGCinfo.txt',header=T,stringsAsFactors=F)
	
	univDmd=univDmd[univDmd$fileId!='',]
	str(univDmd)
	
	getPlot=function(candf){
		candf=candf[candf$relDepth>0,]
		fileIds=levels(factor(candf$fileId))
		
		for(fileId in fileIds){
			fileInfoset=subset(candf,fileId==fileId)
			bmp(paste(fileId,'_readsDepth.bmp',sep=''),height=800,width=1200)
			plot(density(fileInfoset$reads),xlim=c(0,20000),main=fileId)
			dev.off()
			
			}
		
		}
	
	#############################################
	# R scripts ending here
	#############################################
	
	
	now searching on chrX for different GC content area
	note for python filehander, when using seek, the offset starts from 0
	
	
1102:
	!!!!!!!!!!funny fact:
	when using apply_async, pasing args using args=(,), ',' is necessary for becoming a tuple, or a list will be returned
	
	when applying muptiprocessing, variables seems not been shared, this could be solved using Value or Array
	but in-function write file is prefered
	
	get prepared for 11.14 report
	get breakpoint method for 27K, remove the remapping procedure
	get all possible info ready for ML
	get prepared for paper reading
	EWT algos may also work for 27K, adjust for detailed overall error finding rate, this is the only param that matters
	
	ML is just a tool for achieving goal
	
	
	
	34: 249 0825
		/results/analysis/output/Home/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/plugin_out/DMD_plugin_out.1865
		
	35: 615 0831
		/results/analysis/output/Home/Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955/plugin_out/DMD_plugin_old_version_out.6019
	35: 630 0912
		/results/analysis/output/Home/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/plugin_out/DMD_plugin_old_version_out.6021
	35: 695 1014
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/plugin_out/DMD_plugin_out.6009
	35: 697
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/plugin_out/DMD_plugin_out.6040
		
	37: 617 0930
		/results/analysis/output/Home/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/plugin_out/DMD_plugin_out.1594
	37: 621 1008
		/results/analysis/output/Home/Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/plugin_out/DMD_plugin_out.1593
	
	
	! note for python sort, str and int sort has totally different results
	
	
1103:
	get stats out, get all DMD duplication info
	
	for apply_async:
		this would return a exit status which indicate current process processing correctly or not
	
	take time for arranging
	
	
	
1107:
	in formal work, we have calculated large gap statistics according to bed file
	something we should keep in mind:
		for DMD case, duplication and deletion exist for large fragment, thus 200bp bed stats would count
		caming back to our BP project, for heters female, their BP site would have 50% chance to passing to their offsprings, thus testing male deletion would be just enough for calling BP features
	thus for today:
		testing soft-clipping method for all existing file, and then transfer to 45 local directory
		apply original DMD plugin medhod for soft-clipping testing
		thinking about new strategies for data backups
		
	things to do this afternoon:
		1\ get known detailed algos from comparative reviews
		2\ faire des notes sur les methodes et le faire realizer
		
		
	modify the json file and rerun DMD-plugin in 697 run, specially for barcode 79 & 80
	
	
1108:
	get recent DMD data ready
	
	bed file correction:
		
		DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed
		this bed contains both dmd and ref regions, which is already exon annotated
		
	eliminate GC correction procedure for all following work
	
	type I: false positive
	type II:false negative
	
	for EWT:
		false positive ratio is achieved by simulation
		
	one important consideration about chosing window size is that distribution of reads should be approx to be normal distribution
	
	use q20 as a filtration threshold, this would filter out a large part of noises
	
	numpy and scipy do exists in python 2.6.5
	
	#############################
	### EWT PROCEDURE FOR 27K ###
	#############################
	1\ determine different size of window, and plot the reads distribution pattern
	2\ apply no GC adjustment
	3\ assign to little functions
	
1109:
	note that in numpy fancy indexing, the sliced array has always been copied
	
1111:
	finished EWT with numpy
	
	get prepared for DMD paper writing
	
1112:
	a python version of getting returned value from multiprocessing.Pool.apply_async
	
	#########
	import multiprocessing
	import time

	def func(msg):
		print "msg:", msg
		time.sleep(3)
		print "end"
		return "done" + msg

	if __name__ == "__main__":
		pool = multiprocessing.Pool(processes=4)
		result = []
		for i in xrange(3):
			msg = "hello %d" %(i)
			result.append(pool.apply_async(func, (msg, )))
		pool.close()
		pool.join()
		for res in result:
			print ":::", res.get()
		print "Sub-process(es) done."
	##########
	
	watch out for while loop included with a if
	
	
	
1114:
	prepare for DMD paper-writing
	
	use shapiro-test for testing normal distribution: higher P-value would yield a normal distribution pattern
	Kolmogorov-Smirnov (K-S test) were also applied for testing target data with target distribution
	
	for(subfile in dir()){
		data=read.table(subfile,header=F,stringsAsFactors=F)
		data1=data[data[,1]>400,1]
		hist(data1,breaks=100)
		shapiro.test(data1)
		}
		
	use samples with no dup/dels for estimating FPR and window size
	
	first filtra out non-CNV samples
	
	630:
		IonXpress_068
	617:
		IonXpress_016
		IonXpress_006
		IonXpress_007
		IonXpress_011
	621;
		IonXpress_057
		IonXpress_013
		IonXpress_062
		IonXpress_018
		IonXpress_022
	695:
		IonXpress_053
		IonXpress_029
	697:
		IonXpress_057
	646:
		IonXpress_054
		IonXpress_055
		IonXpress_051
		IonXpress_059
		IonXpress_032
		IonXpress_030
		IonXpress_015
	704:
		IonXpress_078
		IonXpress_075
		IonXpress_073
		IonXpress_070
		IonXpress_024
		IonXpress_006
		IonXpress_020
		IonXpress_068
		IonXpress_060
	712:
		IonXpress_053
		IonXpress_007
		IonXpress_023
		IonXpress_063
		IonXpress_062


1115:
	rewrite 27K previous method
	adding 1X 10X 30X 50X coverage info in future report
	PCR dup can not be removed due to depth reduction


	
1116:
	testing 27 method:
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169/dmd_project
		
		
	use np.savetxt('name',obj,delimiter='\t') and np.loadtxt
		instead of using tofile, which would return a flat file
	
	
1117:
	34: 249 0825
		/results/analysis/output/Home/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/plugin_out/DMD_plugin_out.1865
		
		
	35: 615 0831
		/results/analysis/output/Home/Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955/plugin_out/DMD_plugin_old_version_out.6019
	35: 630 0912
		/results/analysis/output/Home/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/plugin_out/DMD_plugin_old_version_out.6021
	35: 695 1014
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/plugin_out/DMD_plugin_out.6009
	35: 697 1026
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/plugin_out/DMD_plugin_out.6040
	35: 704 1104
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-450-P35-HiQ-DMD-32ZA1-Can28-CHMO-CVD-pooling-20161104_704_1055/
	35: 712 1112
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169/
		
		
	37: 617 0930
		/results/analysis/output/Home/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/plugin_out/DMD_plugin_out.1594
	37: 621 1008
		/results/analysis/output/Home/Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/plugin_out/DMD_plugin_out.1593
	37: 646 1028
		/media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646
	
	python cnv_calling_27k_v0.04.py -A /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169 -R /home/ionadmin/ljzhang/DMD_27k_rewrite/712 -U 8
	
	column -t IonXpress_055_rawlib_sorted.bam.dmdRelDpth.reprocessed.results |le
	
	
1118:
	DMD raw data bcp dir:
	DMD raw data backup dir:
	DMD raw data back up dir:
		/media/disk2/ywang/rawdata/rawdata_35/DMD
		/media/disk4/rawdata/rawdata_34/DMD
		/media/disk4/rawdata/rawdata_35/DMD
		/media/disk4/rawdata/rawdata_37/DMD
	
	only need to rewrite bamStatv3.pl using python
	this script contains infos as:
	
		IonXpress_055_rawlib_Basic_statistics_for_mapping_data.txt
	
				target region size bp
				total reads
				total bases
				
				pct of bases with quality >= 20			
				pct of bases with quality >= 30			
				
				mapped reads							using FLAG info
				mapped read rate
				mapped bases							using CIGAR info and reads length
				mapped base rate
				
				on target reads							1bp overlapped reads is treated as target reads
				on target read rate
				on target bases							using CIGAR
				on target base rate
				
				duplicated reads						mapping pos chr:stt:end as key, appearing time as value
				duplicated reads rate					duplicated reads/mapped reads counts
				
				1X coverage of target region			base coverage, using a univ dct
				10X coverage of target region
				20X coverage of target region
				30X coverage of target region
				50X coverage of target region
				100X coverage of target region
				200X coverage of target region
				
				mean depth of target region				calculate with on target base
				
				mismatch reads							original script using optional fields, but CIGAR is prfred
				mismatch read rate
				mismatch base rate
				
				insertion read rate						using CIGAR
				insertion base rate
				
				deletion read rate						using CIGAR
				deletion base rate
				
				eveness score							using mackry et al. 2010 method
				mean read length
				mapped read length
				unmapped read length
				target region read length
	
		
		IonXpress_055_rawlib_coverage_and_mean_depth_for_each_fragment.txt
		
			chr  sttP  endP  1X  10X  20X  30X  meanSeqDpt  relSeqDpt(to whole target)
			
	
	
1121:
	recent results starting from 697, thus rearranging results from that,
	prepare ready pipeline fine for raw-data, open a 3 core process for processing all the necessary files,
		bam stats and other time cosuming processes can be skipped
		
	
	697 646 704 712 report genesis re-run is ready and propompt to local dir
	
	D:\data_chunqiu\20140522:		negative control and patients
	D:\data_chunqiu\20140627:		double blind experiments
	
	
	get re-write report ready before 17:30
	
	one reference for MD tag:
		https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files
	note that MD tag together with CIGAR could get out mutation info out
	
	getting mismatch info out:
		
		my $isMmRead = 0;
			$string =~ s/(?:\^[AGCT]+)//g;#remove deleted bases,and inserted bases are not in the MD tag
			foreach my $i (split //,$string) {#correct on 16th,Dec,2013
				if ($i =~ /[AGCT]/) {
					$isMmRead = 1;
					$mmBaseCnt++;
				}
			}
			if ($isMmRead) {
				$mmReadCnt++;
			}
			
	note that in python, re.sub(regx,tgt,obj) returns str after substitution
		while re.subn() returns a tuple containning after substitution and match times
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	$mappedReadCnt++;
	$mappedBaseCnt += length($SEQ);
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	note that mappedReadCnt is calculated based on FLAG info, which may pose a problem
	
	nX coverage of target region is given as a percentage
	
	
1122:
	a quick tvc pipeline is in running:
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-464-P35-HiQ-DMD-2M-8ZA1-pooling-P17-pooling-TEST-20161119_721_1085/plugin_out/DMD_plugin_out.6215/dmd_project
		/home/ionadmin/ljzhang/project/721_quick_tvc
		
		1123update: also this tvc script uses target bed file which would not apply to 2M experiment
		
	DMD_pluginPipeline is ready on run at:
		/home/ionadmin/ljzhang/20161122.*.out
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-464-P35-HiQ-DMD-2M-8ZA1-pooling-P17-pooling-TEST-20161119_721_1085/plugin_out/DMD_pipline_rerun
		
		1123update: pipeline would not work for this is a 2M experiment, tvc shows no results
	
	
1123:
	721 abd 722 run is online:
	note that 721 is a 2M sample
	
	621 sample also contains 8 2M samples
		use it for analysis
	
	add bamStat function to 2M pipeline
	
	when checking tvc results using AF=AO/DP
		raw vcf would be prior to filtered vcf(which is char by low quality)
		
	
	721 and 722 Zhangcheng DMD sample tvc results are stored at 
		35 server:
		/home/ionadmin/ljzhang/DMD_zhangcheng_single_cell
	
	for 2M DMD breakpoint:
		45 server:
			721 results: /media/disk2/ljzhang/project/dmd/721_20161121_dmd_BP
			722 results: /media/disk2/ljzhang/project/dmd/722_20161121_dmd_BP
			621 results: /media/disk2/ljzhang/project/dmd/621_BP_1008
		
	should merge bamstats function for 2M and 17 pannal
	
	using samtools tview for viewing depth info:
		samtools tview IonXpress_046_rawlib.sorted.bam /results/plugins/IDP3500_plugin/db/db/hg19.fasta
		
	a quick tvc script is ready at 45:
		/home/ionadmin/ljzhang/project/quick_tvc_20161123.py
		
	DMD 2M pipeline is rewrote using pool.apply_async instead of using map
	avoid passing global variables inside a parallel function
	
	
1124:
	USE UNIFIED NAMING METHOD FOR AVOID THE FRAGMENTATION
	get DMD breakpoint samples away from 27K and P16 samples
	get bamStat ready for each pannal, including sample size info
	
	721 722 using M2 labelled and rerun
	
	further alteration is that if it is a partial region capture of whole region capture
	use this for naming standard
	
	
	serveral things not finished:
		1\ CNV calling match
		2\ DMD breakpoint info summary
		3\ bamStats rewrite
		
	note that when we talking about coverage, same thing would apply to depth
	
	721/722 using dmdTest plugin still would not run
	sequencing quality problem
	
	712 provides better samples for testing bamStat:
		
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169/dmd_project
		
		
	python bamStat.py -D /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169/dmd_project -R /home/ionadmin/ljzhang/bamStatRewrite
		
	
20161125:	
	rewrited bamStat and originals differs in on-target and dups
	also target size is greater, thus eliminate ref regions
	note that bamStat.v3.pl uses DMD100_Targets.bed as bed input instead of DMD100_Probes
	
	reClaim that:
		DMD uses transcript of NM_004006
	
	note that in both multi-threading and nulti-processing programing, args is passed down by a tuple, thus ',' is necessary
	
	use 'shutil.rmtree()' for removing trees with files
	
	new DND processing procedure:
		1\ use original z-score method
		2\ consider all other samples in a run as negative samples
		3\ use t-test with negative samples
		4\ use qPCR method, use 2^delta(testing-control)
		5\ for 1-4 steps, use remove duplication filtration
	
20161128:
	pick out a run with no str var or mutations
	
	usng tview for viewing multi-reads alignment to a reference genome
	an example containing complete set of commands are in following:
		$ bwa index reference.fas 
		$ bwa aln reference.fas short_reads.fas >short_reads.sai
		$ bwa samse reference.fas short_reads.sai short_reads.fas >short_reads.sam
		$ samtools faidx reference.fas
		$ samtools import reference.fas.fai short_reads.sam short_reads.bam
		$ samtools sort short_reads.bam short_reads.srt
		$ samtools index short_reads.srt.bam
		$ samtools tview short_reads.srt.bam reference.fas
	
	samtools rmdup -s input.srt.bam out.bam
	
	cnv_calling_27k_v0.05.py uses ref mean number as baseline to calculate relative depth
	
		new DND processing procedure:
		1\ use original z-score method
		2\ consider all other samples in a run as negative samples
		3\ use t-test with negative samples
		4\ use qPCR method, use 2^delta(testing-control)
		
		5\ for 1-4 steps, use remove duplication filtration
		
		
		
	run-selecting:
		695-run:
		35:/results/analysis/output/Home/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037
		35:/media/disk4/rawdata/rawdata_35/DMD/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037
		
		617-run:
		37:/results/analysis/output/Home/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010
		35:/media/disk4/rawdata/rawdata_37/DMD/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010
		
		630-run:
		35:/results/analysis/output/Home/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980
		35:/media/disk4/rawdata/rawdata_35/DMD/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980
		
		249-run:
		34:/results/analysis/output/Home/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376
		35:/media/disk4/rawdata/rawdata_34/DMD/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376
	
	
	DMD raw data bcp dir:
		/media/disk2/ywang/rawdata/rawdata_35/DMD
		/media/disk4/rawdata/rawdata_34/DMD
		/media/disk4/rawdata/rawdata_35/DMD
		/media/disk4/rawdata/rawdata_37/DMD
	
	use samtools view -H IonXpress_018_rawlib.basecaller.bam |head -n 2 |tail -n 1 for getting gender-info
		view SM colomn info, containing sample info
		
	get raw file ready while using R for testing.
		note that large deletion should be removed in order for detecting duplications
		
		? one possible solution is that getting out mean and std and CV out for each run
			after removing large dups & dels, use t test method 
			
	
20161129
	rsync -a ionadmin@192.168.0.35:/results/plugins/DMD_plugin ../data
	
	https://www.ncbi.nlm.nih.gov/genome/guide/human/
	parkinsen
	
	nohup python cnv_calling_27k_test_univ_info_v0.2.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/695 -A /media/disk4/rawdata/rawdata_35/DMD/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/basecaller_results -U 3 >695.log 2>695.nohup.out &

	
	nohup python cnv_calling_27k_test_univ_info_v0.2.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/617 -A /media/disk4/rawdata/rawdata_37/DMD/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/basecaller_results -U 3 >617.log 2>617.nohup.out &
	
	
	nohup python cnv_calling_27k_test_univ_info_v0.2.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/630 -A /media/disk4/rawdata/rawdata_35/DMD/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/basecaller_results -U 3 >630.log 2>630.nohup.out &
	
	
	nohup python cnv_calling_27k_test_univ_info_v0.2_test_error.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/249 -A /media/disk4/rawdata/rawdata_34/DMD/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/basecaller_results -U 3 >249.log 2>249.nohup.out &
	note that in 249 sample head info differs
	
	
20161130
	get all stats file ready
	raw files stored at: D:\thor\documents\WorkProjects\DMD_Project\DMD_pluginRewriteTest\4run_test_data
	
	
	####################### METHODS #########################
	# 1  original z-score
	#   1.1  use zscore method based on inner-sample of different intervals
	#   1.2  use zscore method based on inter-sample
	# 2  applying inter-run zscore
	#   2.1  using intra-sample mean and sd for each interval
	# 3  use t-test based on negative samples
	# 4  use qPCR method, meganify signals using a 2-based exponential distribution
	####################### METHODS #########################
	
	!!!!!!!!!!!!!!!
	we should reconsider the possibility of DMD plugin is wrong, which has applied ref region to a wrong way
	delCoef and dupCoef could be used to align all samples together using Z-score
	!!!!!!!!!!!!!!!
	
	a new processed info head would be:
		
	
	
	
20161201
	paste(vector,sep='',collapse=',')
		http://blog.sina.com.cn/s/blog_6caea8bf0100xiy2.html
	use %in% and setdiff for set computation
	data.frame can be extended
	
	get new DMD plugin out ASAP
	
	
20161102
	use 
	arr[np.ix_([],[])] for subsetting matrix region for analysis
	
	get prepared for paper reading
	
	
	prepare dup removed data for analysis
	nohup python cnv_calling_27k_test_univ_info_v0.2.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/695rd -A /media/disk4/rawdata/rawdata_35/DMD/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/basecaller_results -U 3 -M Y >695rd.log 2>695rd.nohup.out &
	nohup python cnv_calling_27k_test_univ_info_v0.2.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/617rd -A /media/disk4/rawdata/rawdata_37/DMD/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/basecaller_results -U 3 -M Y >617rd.log 2>617rd.nohup.out &
	nohup python cnv_calling_27k_test_univ_info_v0.2.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/630rd -A /media/disk4/rawdata/rawdata_35/DMD/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/basecaller_results -U 3 -M Y >630rd.log 2>630rd.nohup.out &
	nohup python cnv_calling_27k_test_univ_info_v0.2_test_error.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/27K_multiple_test/249rd -A /media/disk4/rawdata/rawdata_34/DMD/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/basecaller_results -U 3 -M Y >249rd.log 2>249rd.nohup.out &
	
	
	python probe100_exon_anno.py /media/disk2/ljzhang/data/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed /media/disk2/ljzhang/project/27K_multiple_test/249rd/249rd.txt
	python probe100_exon_anno.py /media/disk2/ljzhang/data/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed /media/disk2/ljzhang/project/27K_multiple_test/617rd/617rd.txt
	python probe100_exon_anno.py /media/disk2/ljzhang/data/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed /media/disk2/ljzhang/project/27K_multiple_test/630rd/630rd.txt
	python probe100_exon_anno.py /media/disk2/ljzhang/data/DMD_plugin/DMD100_Probes.sorted.merged.4col.based.on.79exon.split.by.200bp.bed /media/disk2/ljzhang/project/27K_multiple_test/695rd/695rd.txt

	
	
	
	for numpy:
		univnpLst[:,2]=univnpLst[:,2].astype(int)/3
		univnpLst is a str np.array
		after division, results are stored back in univnpLst in type of STRING!!
		good news for conducting HPC
		
	
	Reduce(intersect,  list(v1 = c("a","b","c","d"), 
              v2 = c("a","b","e"), 
              v3 = c("a","f","g")))
	for getting intersect between multiple dataset
	
	
1205
	report preparation and folder rearrangement
	
1206
	for 53_20161130 run:
		the dmd breakpoint part:
			ionXpress 12 and 50 has the same point mutation, should exclude the possibility of sample miss-taken
	
	run a tvc pipeline on 12 and 50 to exclude the possibility of sample dups.
	
	a two vcf cmp python script is ready
	note that we should get each script common feature and avoid dup works
	
	
	get a pre-version out with unnecessary options removed
	
	use traceback:
		returns error in a string
			traceback.format_exc()
	
	use fcntl for look file handle
	fcntl only proves 
	
			import fcntl
			f=open('./test')
			fcntl.flock(f,fcntl.LOCK_EX)
			fcntl.flock(f,fcntl.LOCK_UN)
	
	note that fcntl only proves to be efficient under linux which uses a POSIX interface of UNIX
	
	print full path info
	
	
	
1207:
	
	ls IonXpress_*_rawlib.bam | awk '{print "ln -s /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169/dmd_project/",$0}' OFS='' > ../../dmd_plugin_v2_test/dmd_project/2.sh
	
	ls IonXpress_*_rawlib_sorted.bam* |awk '{print "ln -s /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/DMD_plugin_out.6169/dmd_project/",$0}' OFS='' > ../../dmd_plugin_v2_test/dmd_project/3.sh
	
	ls IonXpress_*_rawlib.basecaller.bam | awk '{print "ln -s /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/basecaller_results/",$0}' OFS='' > ../../dmd_plugin_v2_test/dmd_project/1.sh
	
	scp ljzhang@192.168.0.45:/home/ljzhang/project/dmd_plugin_pipeline_v2.0.py ../dmd_plugin_v2_test/
	
	python dmd_plugin_pipeline_v2.0.py -R /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/dmd_plugin_v2_test -B /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/basecaller_results -A /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069 -U 8 -K n
	
	
	python dmd_plugin_pipeline_v2.0.py -R /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/dmd_plugin_v2_test -A /results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069 -M 1 -I IonXpress_003_rawlib.basecaller.bam
	
	
	/results/analysis/output/Home/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/plugin_out/dmd_plugin_v2_test/CNV_calling_tmp/tmpDepthFiles
	
	
1208:
	get erro log of list out of index error of cnv_calling_27k_v1.0.py
	
	scp ljzhang@192.168.0.45:/home/ljzhang/project/dmd_plugin_pipeline_v2.0.py ../beta_dmd_plugin_v2.0/
	scp ljzhang@192.168.0.45:/home/ljzhang/project/cnv_calling_27k_v2.0.py ../beta_dmd_plugin_v2.0/
	scp ljzhang@192.168.0.45:/home/ljzhang/project/launch.sh ../beta_dmd_plugin_v2.0/
	
	
1209:
	write a pipeline for all backed up files:
		change method for getting out gender info, eliminate the pick_up_dmds_procedure
		
	getting prepared for paper reading
	
	
	cp dmd_plugin_pipeline_for_backup_files_v1.0.py /media/disk2/ljzhang/data/DMD_plugin/
	cp cnv_calling_27k_for_backup_files_v1.0.py /media/disk2/ljzhang/data/DMD_plugin/
	cp pick_out_dmd_bams_from_a_run_for_backup_files_v1.0.py /media/disk2/ljzhang/data/DMD_plugin/
	
	testing script:
		python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_for_backup_files_v1.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_test_user_forbcp_test_617 -B /media/disk4/rawdata/rawdata_37/DMD/Auto_test_user_forbcp_test_617/basecaller_results -A /media/disk4/rawdata/rawdata_37/DMD/Auto_test_user_forbcp_test_617 -U 8
	
	
1210
	when appling script to original file, try to make a copy for avoiding rewriting problem, or a entire run data is wasted
	
1212
	avoid erasing data of a backup by making a copy of a chosen example
	get that pipeline of backups going
	
	extract the raw sequence from breakpoint region
	prepare for extracting raw sequence from individual sequencing data
	
	
	
	ll -hl /media/disk2/ywang/rawdata/rawdata_35/DMD/*/basecaller_results/
	ll -hl /media/disk4/rawdata/rawdata_34/DMD/*/basecaller_results/
	ll -hl /media/disk4/rawdata/rawdata_35/DMD/*/basecaller_results/
	ll -hl /media/disk4/rawdata/rawdata_37/DMD/*/basecaller_results/
	
	
	
	note that we ruined two run:
		/media/disk4/rawdata/rawdata_37/DMD/user_sn247560013-44-P37-Hi-Q-T21pooling-Cancer-DMD-20161029_re_1066/basecaller_results/
		
		/media/disk4/rawdata/rawdata_37/DMD/Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/basecaller_results/
		
		
	making a copy of:
		/media/disk4/rawdata/rawdata_37/DMD/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/
		dest:
		/media/disk4/rawdata/rawdata_37/DMD/Auto_user_forbcp_test_617
		
		
		
	cp dmd_plugin_pipeline_for_backup_files_v1.0.py /media/disk2/ljzhang/data/DMD_plugin/
	cp cnv_calling_27k_for_backup_files_v1.0.py /media/disk2/ljzhang/data/DMD_plugin/
	cp pick_out_dmd_bams_from_a_run_for_backup_files_v1.0.py /media/disk2/ljzhang/data/DMD_plugin/
		
		
	python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_for_backup_files_v1.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_forbcp_test_617 -B /media/disk4/rawdata/rawdata_37/DMD/Auto_user_forbcp_test_617/basecaller_results -A /media/disk4/rawdata/rawdata_37/DMD/Auto_user_forbcp_test_617 -U 8
		
		
	
	info summary:
		707 run:
		/results/analysis/output/Home/Auto_sn247560054_sn247560054-451-P35-HiQ-pooling-P17-20161109_707_1057/plugin_out/DMD_plugin_out.6171
		possible info error
		
		712 results ready
		
		713 sequenced as WES
		
		721 run:
			/results/analysis/output/Home/Auto_sn247560054_sn247560054-464-P35-HiQ-DMD-2M-8ZA1-pooling-P17-pooling-TEST-20161119_721_1085/plugin_out/DMD_plugin_out.6215
			
		722 run:
			/results/analysis/output/Home/Auto_sn247560054_sn247560054-465-P35-HiQ-pooling-P17-DMDQ2-20161122_722_1087/plugin_out/DMD_plugin_out.6231
			run failure
			
		728 run:
			/results/analysis/output/Home/Auto_sn247560054_sn247560054-471-P35-HiQ-pooling-P17-dmd19za1-161125_728_1102/plugin_out/DMD_plugin_out.6285/CNV_result
			
		734 run:
			/results/analysis/output/Home/Auto_sn247560054_sn247560054-477-P35-HiQ-P17-DMD-Pooling-21za1-16s-161202_734_1115/plugin_out
			
		53 run: 37 server
			/results/analysis/output/Home/Auto_user_sn247560013-53-P37-HiQ-pooling-dmd-18za1-20161130_53_016/plugin_out/DMD_plugin_out.33
		
		
		previously:
			13 DMD cases with poor quality remains not on machine
			
			
		sq:
			Auto-sn247560054_SN2-432-P35-HiQ-pooling-P17-45ZA1-160923_638_995
			P17_XSQ-31
			
			
	
1213:
	prepare for rerun qPCR for unclear samples
	
	use snownlp for chinese char processing, including pinyin autogenesis
	
	use on-line dataase to see if 
	
	two databases which may be useful for annotating point mutations
		https://www.ncbi.nlm.nih.gov/clinvar/?term=DMD%5Bgene%5D
		https://www.ncbi.nlm.nih.gov/snp/?term=DMD%5Bgene%5D
	
	for taking mutation annotation:
		first check out allele frequency in population(1k freqs)
		if 1k Freq < 0.05:
			check out on clinvar databases
			
			
	when design qPCR analysis, expand one extra exon for both sides:
	for example:
		for exon 45-53 deletion:
			verify exon 44-54
	
1214:	
	try to merge CNV detecting methods
	
1215:
	filter dup and get sample info ready
	filter out a run for pre qPCR CNV verification
	
	gender:
		if hom/total > 0.75 :
			male
		
		
	use column -t file | less -S for align each column
	
	multiple methods for getting out shell pipeline results
	
	
	high_confidence_homozygous_deletion
	high_confidence_homozygous_duplication
	homozygous_duplication
	heterozygous_deletion
	heterozygous_duplication
	
	filehandle can be passed as reference
	
	
1216:
	prepare for all patients summary info
	
	
1219:
	re-talk about the qPCR verifying method
		note that there should be 3 parallels for each sample, thus:
		955 if all exons been verified
		776 if 3 exons span been verified
		533 if 2 exons span been verified
		
	1\ re merge 
	
	
	for unicode problems:
		def test2():
			os.chdir(u'D:\\thor\\documents\\WorkProjects\\DMD_Project\\DMD项目管理\\DMD下机数据汇总报告\\DMD_qPCR后期验证')
			path=os.getcwd()
			path=unicode(path,'gbk')
			print path
			
			lastPat=re.compile(r'.*\\(.+$)')
			lastName=lastPat.findall(path)[-1]
			print 'last name'
			if lastName == u'DMD_qPCR后期验证':
				print lastName
				print 'match'
	
		def test3():
			os.chdir(u'D:\\thor\\documents\\WorkProjects\\DMD_Project\\DMD项目管理\\DMD下机数据汇总报告\\DMD_qPCR后期验证')
			
			infh=open('dmd_nbr_name.txt','r')
			
			names=[]
			univinfo=[]
			for line in infh.xreadlines():
				linear=line.strip().split('\t')
				names.append(unicode(linear[0],'gbk'))
				univinfo.append(line.strip())
				
			for name in names[:10]:
				print name
				
			tgtar=[st for st in univinfo if u'曾灿杰' in unicode(st,'gbk')]
			print unicode(tgtar[0],'gbk')
	
	
	first get a merged script using python, then we seek for best parameter set
		avoid using magic params, 
	
1220:
	get a merged script for CNV calling using new method
	
	get nearby sequence around breakpoint we detect
	
	
1221:
	get homo sequence!!!!
	merge script!!!
	prepare for the annual report
	
	
	caution:
		when writing python, when passing args to functions, do not use mixed indicating method
	
	get out all DMD breakpoint runs
	
	use processed CNV info for analysis
	
	a historical DMD breakpoint project remains unanalysised:
		/media/disk4/rawdata/rawdata_37/DMD/Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/basecaller_results
		not that barcode 3 has mixed info, thus remains unanalysised
		617 is not a breakpoint project file
	
	use Clustal TOOLKIT FOR sequence homology alignment
	
	
1222：
	use refGene (http://varianttools.sourceforge.net/Annotation/RefGene):for different gene reference
	
	use random.sample(list,nbrOfChoice)
	
	for dmd attachment file:
		get depth file ready for each run
		rearrange for experiment
	
	barplot(as.numeric(cnvInfo[,3]),beside=T,names.arg=cnvInfo[,18],cex.names=.6,space=2)
	
	get all info in one table:
		note that table should correspond to original sample sum table info
		use tags for filtration
		
	
1225:
	be careful for 646 run
	
	
	
1227:
a barplot example:
	CCA <-c(3988, 4129, 2409, 7779)
	names(CCA) <- c("ionosphere", "pima", "bupa", "German")
	CCB <-c(3273, 3269, 2318, 5166)
	names(CCB) <- c("ionosphere", "pima", "bupa", "German")
	CC <- cbind(CCA, CCB)
	barplot(t(CC), beside = TRUE,legend = c("Serial", "Parallel"),width = c(100, 100), args.legend = list(x = "topleft", cex=2),ylim = c(0, 8000),cex = 2,cex.axis=2)  
	
	
	
1228:
	use deletion data from 97 patients for plotting dmd exon deletion Pat
	prepare for work-proposal for 2017
	python quick_tvc_v1.0.py -D /home/ionadmin/ljzhang/DMD_zhangcheng_single_cell/750 -P 8
	
	
1229:
	complete the qPCR ggplot program, prepare for getting different controls, including control-m and control-f
	
	when applying AND(&) logic computation, be careful for the order, for the order correspond to tandem application
	
	use read_excel from readxl for read in excel infomation
	
	
0104:
	take a readme file under folders of scripts
	use this json params for tvc calling:
		CHP2.20130920.somatic_lowstringency_pgm_4.0_parameters.json
	
	note that we have take out all base CNV info from bcpFiles
	
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_sn247560054-450-P35-HiQ-DMD-32ZA1-Can28-CHMO-CVD-pooling-20161104_704_1055
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_sn247560054-464-P35-HiQ-DMD-2M-8ZA1-pooling-P17-pooling-TEST-20161119_721_1085
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_sn247560054_sn247560054-465-P35-HiQ-pooling-P17-DMDQ2-20161122_722_1087
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_user_2456267-0357-199-P34-Hi-Q-deaf-pooling-20160909_264_405
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N P34-Hi-Q-dmd34za1-deaf-20160825-Re_386
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_re_982
	python ~/project/cnv_method_testing_data_preparation_20170104.py -N user_sn247560013-44-P37-Hi-Q-T21pooling-Cancer-DMD-20161029_re_1066
	
	
	be careful for text format of uploading file
	
	note that should include a rm-duplicates switch inside the CNV base info extraction
	change script from:
		cnv_calling_27k_test_univ_info_v0.2_test_error.py
	modified script is named as:
		cnv_calling_27k_for_base_info_extraction_with_dup_removed_v0.1.py
		
	
0105:
	get out the qPCR info for the existing files
	merge all cnv testing method together with standard input
	for qPCR qc:
		for sampleSize > 3:
			if delta Ct > 0.25:
				results would not be accepted
				
	qPCR checking:
		161227-DMD-2_data.txt.signalout			checked
		170103-DMD_data.txt.signalout			checked
		20161226-DMD-1_data.txt.signalout		checked
		20161226-DMD-2_data.txt.signalout		checked
		20161227-DMD-1_data.txt.signalout		checked
		qPCR-20161222-DMD-1_data.txt.signalout	checked
		20170107-DMD-3_data.txt.signalout		checked
			DMD_exon51	DM	m	256736.855
		20170107-DMD-1_data.txt.signalout		checked
			bizzard results with many exon dups
		20170107-DMD-2_data.txt.signalout		checked
			DMD_exon46	CC	m	353.155
			DMD_exon46	ZWH	m	95858.424
			DMD_exon46	ZZ	m	90865.595
		170109-DMD_data.txt.signalout			checked
		
		
	paper deadline is late Febrary
	
	for DMD plugin upload format:
		change folder name from barcode to sample ID
		do not use blank space between strings
		change plot type from pdf to png
		
		put all pre-load files inside a ./sampleId/results folder
		
	
	calculate average percent of donner SNP in patients, using file:
		D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\DMD_单细胞\calculate_donner_SNP_proportion_1
		
	
0109:
	GADPH: Chromosome 12: 6,533,927-6,538,374 forward strand. GRCh38:CM000674.2
	which means should share same expression in male and female
	
	merged all signal files together
	
0110:
	for single cell project:
		get univ Sd for patient set
		for each patient:
			get donner SNP distribution
			
	there 49 SNPs left, considered only present in donner
	
0113:
	46,51 probably because of un-specificity of primer
	
	
	https://www.ncbi.nlm.nih.gov/pubmed/?term=876534&report=abstract&format=text
	use this for downloading paper info
	
	use clinvar and annovar for annotation and prediction of pathogenity
	
	ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar
	hash info for getting paper ID and paper 
	
	variation clinical significance: first clinvar, then freq(threshold set to 0.05), then protein annotation, >2 patho) otherwise set to unsignificant
	
	
	############################# important criterion for annotating pathogenic mutations #############################
	for detail:
		if mixed: non-significant
		if likely mixed: non-significant
		if patho is mixed: pathogenic
		if likely present: then is likely to be pathogenic
		
	hash table is likely stored in:
		tab_delimited/variantion_sum.txt 
		should take further sight in
	
	############################# important criterion for annotating pathogenic mutations #############################
	
0116:
	prepare for the final report of qPCR result
	get ready for CNV paper writing
	get know that annovar.pl script
	
	all qPCR results are ready, thus getting final report ready for update
	
	
	
	Traceback (most recent call last):
	  File "/results/plugins/dmdTestDonnotRunMe/cnv_calling_27k_v2.0.py", line 597, in <module>
		main()
	  File "/results/plugins/dmdTestDonnotRunMe/cnv_calling_27k_v2.0.py", line 134, in main
		os.chdir(analyDir)
	OSError: [Errno 2] No such file or directory: '/results/analysis/output/Home/Auto_sn247560054_sn247560054-516-P35-HiQ-P17-3za1-WES1za1-Pooling-170116_778_1196/dmd_project'
	python /results/plugins/dmdTestDonnotRunMe/cnv_calling_27k_v2.0.py -A /results/analysis/output/Home/Auto_sn247560054_sn247560054-516-P35-HiQ-P17-3za1-WES1za1-Pooling-170116_778_1196 -R /results/analysis/output/Home/Auto_sn247560054_sn247560054-516-P35-HiQ-P17-3za1-WES1za1-Pooling-170116_778_1196/plugin_out/dmdTestDonnotRunMe_out.6597 -U 8
			starts...
	Traceback (most recent call last):
	  File "/results/plugins/dmdTestDonnotRunMe/dmd_plugin_pipeline_v2.py", line 619, in <module>
		main()
	  File "/results/plugins/dmdTestDonnotRunMe/dmd_plugin_pipeline_v2.py", line 178, in main
		os.chdir(resultsDir+'/tmpDepthFiles')
	OSError: [Errno 2] No such file or directory: '/results/analysis/output/Home/Auto_sn247560054_sn247560054-516-P35-HiQ-P17-3za1-WES1za1-Pooling-170116_778_1196/plugin_out/dmdTestDonnotRunMe_out.6597/tmpDepthFiles'
	
0117:
	DMD raw data bcp dir:
	DMD raw data backup dir:
	DMD raw data back up dir:
		/media/disk2/ywang/rawdata/rawdata_35/DMD
		/media/disk4/rawdata/rawdata_34/DMD
		/media/disk4/rawdata/rawdata_35/DMD
		/media/disk4/rawdata/rawdata_37/DMD
	
	?
	what kd and other stand for
	
	ll /media/disk2/ywang/rawdata/rawdata_35/other
	ll /media/disk2/ywang/rawdata/rawdata_37/other
	ll /media/disk2/ywang/kd_result ?
	
	ll /media/disk4/rawdata/rawdata_34/DMD
	ll /media/disk4/rawdata/rawdata_35/DMD
	ll /media/disk4/rawdata/rawdata_37/DMD
	
#	this following script returns all folders under each backup folder
	ll /media/disk2/ywang/rawdata/rawdata_35/other/
		Auto_sn247560054_sn247560054-496-P35-HiQ-Pooling-P35-20161222-cancer_758_1156/
		Auto_sn247560054_sn247560054-498-P35-HiQ-T21-Pooling-20161226-1_760_1160/
		Auto_sn247560054_sn247560054-499-P35-HiQ-Pooling-Can28-deaf-NBS-20161227_761_1162/
		Auto_sn247560054_sn247560054-502-P35-HiQ-Deaf-Pooling-CHMO_2017-01-04_764_1168/
		Auto_sn247560054_sn247560054-503-P35-HiQ-T21-Pooling-20170105-1_765_1170/
	
	ll /media/disk2/ywang/rawdata/rawdata_37/other
		Auto_user_sn247560013-68-P37-HiQ-T21-Pooling-20161221-2_69_049/
		Auto_user_sn247560013-70-P37-HiQ-T21-Pooling-20161229-2_72_053/
		Auto_user_sn247560013-71-P37-HiQ-T21-Pooling-20161229-1_71_055/
		
	ll /media/disk4/rawdata/rawdata_34/DMD
		Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/
		Auto_user_2456267-0357-199-P34-Hi-Q-deaf-pooling-20160909_264_405/
		P34-Hi-Q-dmd34za1-deaf-20160825-Re_386/
		
	ll /media/disk4/rawdata/rawdata_34/IDP17
		Auto_user_2456267-0357-162-P34-Hi-Q-P17-38ZA1_pooling_20160728_222_325
		Auto_user_2456267-0357-180-P34-Hi-Q-P17-PKU-pooling-160822_243_361
		
	ll /media/disk4/rawdata/rawdata_35/DMD
		Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955/
		Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/
		Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/
		sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_re_982/
		
	ll /media/disk4/rawdata/rawdata_35/IDP17
		Auto_sn247560054_SN2-432-P35-HiQ-pooling-P17-45ZA1-160923_638_995
		
	ll /media/disk4/rawdata/rawdata_37/DMD
		Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/
		Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/
		user_sn247560013-44-P37-Hi-Q-T21pooling-Cancer-DMD-20161029_re_1066/
		
	ll /media/disk4/rawdata/rawdata_37/IDP17
		Auto_user_sn247560013-38-P37-HiQ-pooling-P17-40ZA1-20161020_640_1052
		
		
	
	
	a 249 re-run script is as follows:
		python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_for_backup_files_v1.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_test_user_forbcp_test_617 -B /media/disk4/rawdata/rawdata_37/DMD/Auto_test_user_forbcp_test_617/basecaller_results -A /media/disk4/rawdata/rawdata_37/DMD/Auto_test_user_forbcp_test_617 -U 8
	
	
	249 left only:
		/media/disk2/ljzhang/project/27K_multiple_test/249
		for:
			/media/disk4/rawdata/rawdata_34/DMD/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376
	
	
	
0119:
	cut -f 3,4,5,18,19 ./IonXpress_028/CNV_results/IonXpress_028.cnv.results | column -t | less -S
	grep 'pathogenic' ./IonXpress_028/tvc_out/IonXpress_028_annotation.igmd | less -S
	
	ll /media/disk2/ywang/rawdata/rawdata_35/other/
		Auto_sn247560054_sn247560054-496-P35-HiQ-Pooling-P35-20161222-cancer_758_1156/
		Auto_sn247560054_sn247560054-498-P35-HiQ-T21-Pooling-20161226-1_760_1160/
		Auto_sn247560054_sn247560054-499-P35-HiQ-Pooling-Can28-deaf-NBS-20161227_761_1162/
		Auto_sn247560054_sn247560054-502-P35-HiQ-Deaf-Pooling-CHMO_2017-01-04_764_1168/
		Auto_sn247560054_sn247560054-503-P35-HiQ-T21-Pooling-20170105-1_765_1170/
		
	ll /media/disk2/ywang/rawdata/rawdata_37/other
		Auto_user_sn247560013-68-P37-HiQ-T21-Pooling-20161221-2_69_049/
		Auto_user_sn247560013-70-P37-HiQ-T21-Pooling-20161229-2_72_053/
		Auto_user_sn247560013-71-P37-HiQ-T21-Pooling-20161229-1_71_055/
		
	ll /media/disk4/rawdata/rawdata_34/DMD
		Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/
		Auto_user_2456267-0357-199-P34-Hi-Q-deaf-pooling-20160909_264_405/
		P34-Hi-Q-dmd34za1-deaf-20160825-Re_386/
		
	ll /media/disk4/rawdata/rawdata_35/DMD
		Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955/
		Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/
		Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/
		sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_re_982/
		
	ll /media/disk4/rawdata/rawdata_37/DMD
		Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/
		Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/
		user_sn247560013-44-P37-Hi-Q-T21pooling-Cancer-DMD-20161029_re_1066/
	
	
0120:
	get dir structure for dmd plugin pipeline script and its subscript
	try to unify script for plugin and for historical backup files
	try to move data to a upload folder
	prepare for point mutation clinical significance prediction
	
	next:
		annova.pl for different clinical prediction
	
	46,51 probably because of un-specificity of primer
	
	https://www.ncbi.nlm.nih.gov/pubmed/?term=876534&report=abstract&format=text
	use this for downloading paper info
	
	use clinvar and annovar for annotation and prediction of pathogenity
	
	ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar
	hash info for getting paper ID and paper 
	
	variation clinical significance: first clinvar, then freq(threshold set to 0.05), then protein annotation, >2 patho) otherwise set to unsignificant
	
	############################# important criterion for annotating pathogenic mutations #############################
	for detail:
		if mixed: non-significant
		if likely mixed: non-significant
		if patho is mixed: pathogenic
		if likely present: then is likely to be pathogenic
		
	hash table is likely stored in:
		tab_delimited/variantion_sum.txt
		should take further sight in
	############################# important criterion for annotating pathogenic mutations #############################
	
	ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt
	
	this link stores all allele ID together with its citation_ID
		AlleleID	VariationID	rs	nsv	citation_source	citation_id
	rs ID would be important for saving from constructing an additional hash table
	
	use this for getting out paper 
	https://www.ncbi.nlm.nih.gov/pubmed/?term=12030328&report=docsum&format=text
	
	the paper annotation woud contain two files:
		1\ rs ID linked to a paper info string
		2\ alleleID linked to a paper info string
	
	
0121:
	133 papers were extracted from ncbi related to DMD & BMD
	
	NOW:
		connect these paper to rs IDs and for preparation of mutation annotation
		prepare for CNV calling method and DMD plugin for existing bams, using relative fewer bam run
		
	
0122:
	get CNV running with V3.0 version
	
	python ~/project/cnv_calling_27k_v3.0.py -P /media/disk2/ljzhang/data/DMD_plugin -U 1 -R /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/20170121_dmdPlugin_test_out/plugin_out/DMD_test.out
	
	note that decorator would not be called inside a multiprocessing module
	
	cp ~/project/cnv_calling_27k_v3.0.py /media/disk2/ljzhang/data/DMD_plugin
	cp ~/project/dmd_plugin_pipeline_v3.0.py /media/disk2/ljzhang/data/DMD_plugin	
	cp ~/project/pick_out_dmd_bams_from_a_run.py /media/disk2/ljzhang/data/DMD_plugin
	
	python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_v3.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/20170121_dmdPlugin_test_out/plugin_out/DMD_test.out -B /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/20170121_basecaller_test_for_dmdPlugin/basecaller_results -A /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/20170121_dmdPlugin_test_out -U 8
	
	
0204
	several things should be accomplished:
		get exon dup/del info for 2M sequencing
		get R plotting script out for hist plot and coverage line
		
	soft-clipping adjusting method should be put afterwards
	
	a summary of all available data should be ready, plus formal backupped file stored at:
		/media/disk2/ljzhang/37_backups/ljzhang (containning complete plotting scripts):
			617_1594out
			621_1593out
			re_20161029_646
		
	
	use old plugin, rearrange results in a single folder
	
	#### paper preparation DMD breakpoint
	D:\thor\documents\WorkProjects\DMD_Project\DMD_breakpoint\2M断点查找项目
		containning run info for 53 and 617
		53 has 2 biased sample while whole 617 used wrong probe(results presented in overall summary)
		
	tips:
		a=sample(1:10,100,replace=T)
		as.data.frame(a)
		aggregate(a,list(a),length)
	
	note that BamStat.v3.a.2015.5.29.pl has generated script used for plotting
		/media/disk2/ljzhang/37_backups/ljzhang/617_1594out/bam_stat_second
		this plot has Mean_depth Coverage and GC_content
		this would be put into front-end presentation
		
	one full 27K results remains in: /media/disk2/ljzhang/37_backups/ljzhang/617_1594out/dmd_project
	
	running to test if any procedure in R plotting is wrong:
		nohup python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_v3.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out -A /media/disk2/ljzhang/project/dmd_617_27k_test -U 4 >0204.log 2>0204.nohup.log &
	
	future experiment plan is stored at:
		D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\2017_02_dmd_paper_preparation
	
	for info crawler part:
		var_citation.txt:
			alleleID	citationID	rsID
		submission_summary_dmd.txt
			variationID	clinicalSignificance
		rs_paper_annotation_dmd.txt
			rsID	citationID	paperRef
		rsID_citation_dict.txt
			rsID	paperRef
		
		use this infomation for point mutation annotation
		
	
0205:
	note that in original python pipeline, R plotting can not be performed normaly, thus should start a single test.
	mainly problem is:
		convert2annovar.pl usually takes too long
		BamStat.v3.a.2015.5.29.pl would not plot the png as we expected
	
	
	objectif:
		finish the function test and fix the timming and plotting problem
		prepare for clusto and blast toolkit
		
	
	use IonXpress_003_rawlib_sorted.bam 			(from 617 test run)
	running command is:
		perl /media/disk2/ljzhang/data/DMD_plugin/BamStat.v3.a.2015.5.29.pl -regionFile /media/disk2/ljzhang/data/DMD_plugin/DMD.enst7033.exon.bed -bamFile IonXpress_003_rawlib_sorted.bam -outDir /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003 -plot -gcCorrection
	
	
	for plotting:
		in original plugin, plugin were assigned to two pdf format:
			IonXpress_003_rawlib_sorted_depthPlot.pdf
			IonXpress_003_rawlib_sorted_depthPlot_with_GC_biase_correction.pdf
			if Cairo package is not installed, then use a pdf_png convertor
			
		gc correction seems to be relatively good, thus adding GC correction in final CNV testing
		seems things would not work with while simply change the suffix
		
		for better presentation:
			add a png function in BamStat.v3.a.2015.5.29.pl for both with and without gc correction
			
	softwarePack using ./convert to convert pdf file to png
		http://www.imagemagick.org
		which presented in original plugin is:
			system("$convert -density 64 $outDir/${sampleName}_depthPlot.pdf $outDir/${sampleName}_depthPlot.png");
			
		exportCmd='export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s'	%	(pluginPathMe+'/tvc_lib')
	
	################### this part aims at solving timming problems ###################
	####
	IonXpress_003_rawlib_sorted.bam
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/media/disk2/ljzhang/data/DMD_plugin/tvc_lib
	/media/disk2/ljzhang/data/DMD_plugin/tvc --parameters-file /media/disk2/ljzhang/data/DMD_plugin/variant_caller_scripts/json/startplugin.mlld2.from.second.tvc4.for.hom.json --reference /media/disk2/ljzhang/data/DMD_plugin/hg19/hg19.fasta --input-bam /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/dmd_project/IonXpress_003_rawlib_sorted.bam --target-file /media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.extend.200bp.flanking.chr.bed --output-dir /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1 --output-vcf IonXpress_003_rawlib.vcf
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib.vcf > /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib.avinput
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/table_annovar.pl -outfile /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib -protocol refGene,avsnp142,1000g2014oct_all,1000g2014oct_eas,1000g2014oct_sas,esp6500siv2_all,exac03,clinvar_20140929 -buildver hg19 -operation g,f,f,f,f,f,f,f -nastring - -remove -otherinfo /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib.avinput /media/disk2/ljzhang/data/DMD_plugin/annovar/humandb/
	perl /media/disk2/ljzhang/data/DMD_plugin/format_annovar_anno.pl /media/disk2/ljzhang/data/DMD_plugin/hg19_Homo_sapiens.GRCh37.63.DMD.enst7033.cds.with.cDNA.numbering.bed /media/disk2/ljzhang/data/DMD_plugin/hg19/hg19.fasta /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib.hg19_multianno.txt /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib_annotation.txt
	perl /media/disk2/ljzhang/data/DMD_plugin/DMD_inhouse_annotate.20150602.pl /media/disk2/ljzhang/data/DMD_plugin/igmd.omim.DMD.txt /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib_annotation.txt /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib_annotation.igmd
	perl /media/disk2/ljzhang/data/DMD_plugin/hgmd.dmd.annotate.pl /media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_substitution.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_splicing.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_small_deletion.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_small_insertion.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_small_indel.txt /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib_annotation.igmd /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out_1/IonXpress_003_rawlib_annotation.igmd.hgmd
	#########################################
	this works fast and solid for 003, also for others
	
	a refined version is running with 8 cores
		(version changed from v3.0 to v4.0)
	
		nohup python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_v4.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out -A /media/disk2/ljzhang/project/dmd_617_27k_test -U 8 >0204.log 2>0204.nohup.log &
	
	
0206:
	note that the time consuming procedure still exits
	these five barcodes were:
		006
		013
		019
		003
		005
		this has nothing related to file size
		and chunk in convert2annovar process
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_006/tvc_out/IonXpress_006_rawlib.vcf > /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_006/tvc_out/IonXpress_006_rawlib.avinput
	
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_013/tvc_out/IonXpress_013_rawlib.vcf > /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_013/tvc_out/IonXpress_013_rawlib.avinput
	
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_019/tvc_out/IonXpress_019_rawlib.vcf > /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_019/tvc_out/IonXpress_019_rawlib.avinput
	
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out/IonXpress_003_rawlib.vcf > /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_003/tvc_out/IonXpress_003_rawlib.avinput
	
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_005/tvc_out/IonXpress_005_rawlib.vcf > /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out/IonXpress_005/tvc_out/IonXpress_005_rawlib.avinput
	
	############
	solution:
		put all annovar script in a single shell and excute
		
	objectif:
		change dmd_plugin_pipeline_v4.0.py for the annovar part
		getout file format for cnv_calling_27k_v3.0.py
		change breakpoint script for more soft-clipping reads
		
	cd /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out
		nohup python /media/disk2/ljzhang/data/DMD_plugin/dmd_plugin_pipeline_v4.0.py -P /media/disk2/ljzhang/data/DMD_plugin -R /media/disk2/ljzhang/project/dmd_617_27k_test/plugin_out/DMD_plugin_out -A /media/disk2/ljzhang/project/dmd_617_27k_test -U 8 >0204.log 2>0204.nohup.log &
	
	
	convert2annovar still takss more time than expected, maybe should keep annovar and bamStat away from multi-core
	
0207:
	strategie change:
		use old plugin while introduce in new CNV calling method
		because annovar takes too long with no reason
			this should be a problem
	
	objectif:
		getout file format for cnv_calling_27k_v3.0.py
		change breakpoint script for more soft-clipping reads
	
	convert2annovar.pl in annovar toolkit may be too old and have bugs,
	
0208:
	objectif:
		get mechanism paper for different gene rearrangement:
		prepare a aim for dmd paper
			two strategies:
				1> get 27K results ready and prepare for more figures and forms, mainly focus on cnv calling method
				2> use dmd 2M results for preparation of breakpoint mechanism analysis
				
0223:
	finished all DMD paper, along with info sorting and categories
	
0302:
	 a sample for patent:
		http://www2.soopat.com/Patent/201410048518
	
	objectif:
		1\ abstract in chinese
		2\ adding breakpoints into methods
		3\ arrange sanger validation
	
	use primer3 for primer design:
		192.168.0.47
			/home/ljzhang/project/dmd/sangerVerif/shell.sh
		use primer premier for validation
	
	
0303:
	serveral modification about the DMD paper:
		1\ correlation of intron length with the enrichment of breakpoints
		2\ reconstruction of all potential LCRs within all DMD intron
		3\ the density of LCR with its effects to breakpoint enrichments
		4\ syntex error & capital error
		5\ adding p value to any conclusion
		4\ arrange 2M experiments for DMD paper, randomize the breakpoint choice.
	
	
0306:
	1\ getting out breakpoint varification table
	2\ dicuss and change if necessary, the landscape of paper
	3\ adding schema of breakpoint theory
	
	
	kind of emergency:
		DMD results has been returned
		results should be cross-verified
		1\ check all results and remap with raw results paper
		2\ get raw info ready for check
		3\ arrange 2M samples according to the new final results
		
		mainly range for deletion
		save duplication to last
		
	manage time well
	
	working flow:
		1\ compare original data with feedbacks
			D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\27K上机汇总信息
				博奥木华DMD基因检测结果核对20170305--中山一DMD结果.xlsx
				dmd_中山医科研样本结果汇总_0118_2017_博奥木华.xlsx
		
	
	rearrangement for file position:
		/media/disk2/ljzhang/37_backups/ljzhang
			617
			621
			re_20161029_646
			
		/media/disk2/ljzhang/project/dmd_617_27k_test
			617 dmd_617_27k_test
			
		/media/disk2/ljzhang/project/27K_multiple_test/
			this is mainly depth info
			249
			617
			630
			695
			
		/media/disk2/ljzhang/project/dmd
			0707_bam_test
			0829_bam_test
			20160831_615_955
			721_20161121_dmd__all
			617	./DMD_cnv_calling/37_617/
			
		/media/disk2/ljzhang/project/DMD_plugin/DMD_27k_algos/34_targetedSam/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/
			249
		
		these 3 were for depth info summary
		/media/disk2/ljzhang/project/DMD_plugin/DMD_27k_algos/35_targetedSam/Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955
			615
		/media/disk2/ljzhang/project/DMD_plugin/DMD_27k_algos/35_targetedSam/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037
			695
		/media/disk2/ljzhang/project/DMD_plugin/DMD_27k_algos/35_targetedSam/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980
			630
		
		/media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin
			Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955
				615
			Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980
				630
			Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037
				695
			Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041
				697
			Auto_sn247560054_sn247560054-450-P35-HiQ-DMD-32ZA1-Can28-CHMO-CVD-pooling-20161104_704_1055
				704
			Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069
				712
			Auto_sn247560054_sn247560054-464-P35-HiQ-DMD-2M-8ZA1-pooling-P17-pooling-TEST-20161119_721_1085
				721
			Auto_sn247560054_sn247560054-465-P35-HiQ-pooling-P17-DMDQ2-20161122_722_1087
				722
			Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376
				249
			Auto_user_2456267-0357-199-P34-Hi-Q-deaf-pooling-20160909_264_405
				264
			Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010
				617
			Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018
				621
			P34-Hi-Q-dmd34za1-deaf-20160825-Re_386
				386
			sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_re_982
				982
			user_sn247560013-44-P37-Hi-Q-T21pooling-Cancer-DMD-20161029_re_1066
				1066
		
		D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\1028上交初汇总\dmd深度图重检验
			621
			646
			695
			697
			704
			
		D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\DMD_qPCR后期验证\dmd_备份文件DMD重跑结果
		D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\tmp
			30
			31
			
		with /media/disk4 data excluded
		
		
	##################
	D:\thor\documents\WorkProjects\DMD文章项目\2017_02_dmd_paper_preparation
		0306结合中山医结果整理报告精简.csv
	

0307:
	1\ getting out 10 cases and cross validate
	2\ arrange for sanger primer
	
	ll /media/disk4/rawdata/rawdata_34/DMD
		Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/						f				
		Auto_user_2456267-0357-199-P34-Hi-Q-deaf-pooling-20160909_264_405/						g
		P34-Hi-Q-dmd34za1-deaf-20160825-Re_386/													g
		
	ll /media/disk4/rawdata/rawdata_35/DMD
		Auto_sn247560054_SN2-413-P35-Hi-Q-DMD-16ZA1-pooling-KD-pooling-20160831_615_955/		g
		Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/							g, sorted included
		Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/					g, sorted included
		sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_re_982/								g
		
	ll /media/disk4/rawdata/rawdata_37/DMD
		Auto_user_sn247560013-18-P37-HIQ-pooling-WES-1-DMD-16ZA1-20160930_617_1010/				g
		Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/		4,6,18,22,50,51,64,66 f
		user_sn247560013-44-P37-Hi-Q-T21pooling-Cancer-DMD-20161029_re_1066/					f
		
	
	use 630 and 695 for report
	
	
	
	samtools sort IonXpress_075_rawlib.bam -o IonXpress_075_rawlib_sorted.bam
	samtools index IonXpress_075_rawlib_sorted.bam
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/disk2/ljzhang/data/DMD_plugin/tvc_lib
	/data/disk2/ljzhang/data/DMD_plugin/tvc --parameters-file /data/disk2/ljzhang/data/DMD_plugin/variant_caller_scripts/json/startplugin.mlld2.from.second.tvc4.for.hom.json --reference /data/disk2/ljzhang/data/DMD_plugin/hg19/hg19.fasta --input-bam IonXpress_075_rawlib_sorted.bam --output-dir /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/IonXpress_075 --output-vcf IonXpress_075_rawlib.vcf

	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/convert2annovar.pl -format vcf4 IonXpress_075_rawlib.vcf > IonXpress_075_rawlib.avinput
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/table_annovar.pl -outfile IonXpress_075_rawlib -protocol refGene,avsnp142,1000g2014oct_all,1000g2014oct_eas,1000g2014oct_sas,esp6500siv2_all,exac03,clinvar_20140929 -buildver hg19 -operation g,f,f,f,f,f,f,f -nastring - -remove -otherinfo IonXpress_075_rawlib.avinput /media/disk2/ljzhang/data/DMD_plugin/annovar/humandb/
	perl /media/disk2/ljzhang/data/DMD_plugin/format_annovar_anno.pl /media/disk2/ljzhang/data/DMD_plugin/hg19_Homo_sapiens.GRCh37.63.DMD.enst7033.cds.with.cDNA.numbering.bed /media/disk2/ljzhang/data/hg19/ucsc.hg19.fasta IonXpress_075_rawlib.hg19_multianno.txt IonXpress_075_rawlib_annotation.txt
	perl /media/disk2/ljzhang/data/DMD_plugin/DMD_inhouse_annotate.20150602.pl /media/disk2/ljzhang/data/DMD_plugin/igmd.omim.DMD.txt IonXpress_075_rawlib_annotation.txt IonXpress_075_rawlib_annotation.igmd
	perl /media/disk2/ljzhang/data/DMD_plugin/hgmd.dmd.annotate.pl /media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_substitution.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_splicing.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_small_deletion.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_small_insertion.txt:/media/disk2/ljzhang/data/DMD_plugin/HGMD_profession_processed_small_indel.txt IonXpress_075_rawlib_annotation.igmd IonXpress_075_rawlib_annotation.igmd.hgmd
	

	cut -f 1-8,18 /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/IonXpress_001/CNV_results/IonXpress_001.cnv.results | less -S

	
	
	
	
0309:
	objectif:
	1\ range aqua folder
	2\ get out diff info betw origin and feedbacks
	3\ prepare for dataStruct and algo writing
	4\ prepare for swiss-waterman and blast algos
	
	
	!!## get out info 
	other DMD bam backup folder: /media/disk5/rawdata/rawdata_35/DMD
	
	grep 'DMD' IonXpress_004/tvc_out/IonXpress_004_annotation.igmd | less -S
	cut -f 1-8,18 ./IonXpress_066/CNV_results/IonXpress_066.cnv.results | less -S

	646 still needs a rerun for tvc and CNV
	
	30 mismatch, in which 14 were confirmed correct, 16 wrong
	163 match and mismatch (133 matched)
	thus correct ratio is (133+14)/163
	
	
0313:
	concentrate on paper work
	objectif of this week
	1\ get out QC script for DMD 2M pipeline
	2\ modify breakpoint caller by removing additional error log, use redirect method and traceback.print_exc()
	3\ get paper ready based on discussion results
	
	objectif today:
		for script modification:
			changed script stored as:
				dmd_breakpoint_caller_v2.21.py
				note that dmd_breakpoint_caller_v2.20.py is still in use for error log option is not in use for original script
		
		for QC:
			27K bamStat is stored as bamStat.py
			2M bamstat named as bamStat_2M.py
			
			bamStata.py function:
			1\ use mapped bam as input, sort would be performed is _sorted.bam do not exist
			2\ use exon table to form a bedDct, should refer to 2M bed file
			3\ change sorted bam to sam file
			4\ qc item includes:
				total reads
				total bases
				
				base > 20
				base > 30
				
				mapped reads
				mapped bases
				
				ontarget reads
				ontarget bases
				
				duplicate reads
				
				x1 bases
				x10 bases
				x20 bases
				x30 bases
				x50 bases
				x100 bases
				x200 bases
				
				mis-match reads
				mis-match bases
				
				insert reads
				insert bases
				
				delete reads
				delete bases
				
				mapped reads length
				mapped reads total length
				
				unmapped read length
				unmapped read total length
				unmapped reads
				
				target read length
				
			note that we could modify bed file to get a fast script of 2M bamStat_2M.py
			
0314:
	a new 2M results is ready on 36 server:
	Auto_sn247770271_DX--160-P36-Hi-Q-DMD-21ZA1-170313_279
		D:\thor\documents\WorkProjects\DMD_Project\DMD项目管理\DMD下机数据汇总报告\2M断点上机汇总信息\170313-ZSY-DMD-21ZA1-上机信息.xlsx
	
	pick out dmd bams from 36 server:
		/home/ionadmin/ljzhang/project/dmd/pick_dmd_samples.py
	
0315:
	recalibrate data
	prepare for report of DMD paper
	contact with Xinzhi for arrangement of DMD experiment, for rearrangement of diff results and for not arranged ones
	make plan for financial reports, < 24w available
	range for materials
	
0316:
	1> get out finacial budget for DMD association
	2> get out CNV and z-score for bioinformatic problematics
	3> get out point mutations for all samples without remarks for dup or dels

	a. get out point mutations
	
	igmd header for igmd file:
		D:\thor\documents\WorkProjects\igmd.head
		
	# this following commands are for getting out DMD point mutations
	# note that a sum-up file is ready and recalibrated in DMD paper project
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/IonXpress_051/tvc_out/IonXpress_051_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/IonXpress_075/tvc_out/retestannovar/IonXpress_075_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/IonXpress_055/tvc_out/IonXpress_055_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/IonXpress_058/tvc_out/IonXpress_058_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_2456267-0357-185-P34-Hi-Q-dmd34za1-deaf-20160825_249_376/IonXpress_053/tvc_out/IonXpress_053_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_user_sn247560013-21-P37-Hi-Q-pooling-DMD-14za1-CVD-CANCER-20161008_621_1018/IonXpress_062/tvc_out/IonXpress_062_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_SN2-425-P35-HiQ-DEAF-Pooling-160909R_630_980/IonXpress_067/tvc_out/IonXpress_067_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_015/IonXpress_015_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_019/IonXpress_019_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_021/IonXpress_021_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_022/IonXpress_022_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_030/IonXpress_030_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_049/IonXpress_049_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/37_backups/ljzhang/re_20161029_646/DMD_plugin_out.1654/dmd_project/quickAnno/tvcAnno/IonXpress_052/IonXpress_052_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-441-P35-HiQ-DEAF-Pooling-161014_695_1037/IonXpress_021/tvc_out/IonXpress_021_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/IonXpress_012/tvc_out/IonXpress_012_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/IonXpress_057/tvc_out/IonXpress_057_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/IonXpress_062/tvc_out/IonXpress_062_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/IonXpress_004/tvc_out/IonXpress_004_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-443-P35-HiQ-DMD-32ZA1-pooling-cancer-CHMO-PKU-20161026_697_1041/IonXpress_005/tvc_out/IonXpress_005_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-450-P35-HiQ-DMD-32ZA1-Can28-CHMO-CVD-pooling-20161104_704_1055/IonXpress_004/tvc_out/IonXpress_004_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-450-P35-HiQ-DMD-32ZA1-Can28-CHMO-CVD-pooling-20161104_704_1055/IonXpress_017/tvc_out/IonXpress_017_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-450-P35-HiQ-DMD-32ZA1-Can28-CHMO-CVD-pooling-20161104_704_1055/IonXpress_028/tvc_out/IonXpress_028_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/IonXpress_063/tvc_out/IonXpress_063_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/IonXpress_007/tvc_out/IonXpress_007_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
		grep 'DMD' /media/disk2/ljzhang/project/dmd_bcpFiles_rerun_with_plugin/Auto_sn247560054_sn247560054-456-P35-HiQ-DMD-34ZA1-pooling-CHMO-can28-20161112_712_1069/IonXpress_071/tvc_out/IonXpress_071_annotation.igmd | cut -f 2,9,14,15,16,17,18-21 | less -S
		
	# point mutations ends here
		
0323:
	table_annovar.pl --nastring - --outfile ./test.anno --remove --otherinfo --buildver hg19 --protocol refGene,snp137,1000g2014oct_all,1000g2012apr_asn,esp6500si_all,exac03,ljb26_all --operation g,f,f,f,f,f,f ./test.anno /home/ybliu/projects/CerebralVascular/bin_v0.2/data/humandb
	
	note that *ljb* has been transfered into /media/disk2/ljzhang/data/DMD_plugin/annovar/humandb
	
	note that table_annovar.pl from DMD_plugin should be modified, for this does not recognize *ljb* reportory
	
	# for sift and polyphen annotation, use command below
	perl /media/disk2/ljzhang/data/DMD_plugin/annovar/table_annovar.pl -outfile IonXpress_017_rawlib -protocol refGene,avsnp142,1000g2014oct_all,1000g2014oct_eas,1000g2014oct_sas,esp6500siv2_all,exac03,clinvar_20140929,ljb26_all -buildver hg19 -operation g,f,f,f,f,f,f,f,f -nastring - -remove -otherinfo IonXpress_017_rawlib.avinput /media/disk2/ljzhang/data/DMD_plugin/annovar/humandb/
	
	
0324
	get out budget for sanger and 27K and send out for validation
	valuation would be carried out next week
	should note in the budget, only 80% percent is now available, budget should be made till augest
	
0327:
	note that in budget, there remains 184 k for testing dosage
	split the budget to 4 months
	
	review about LCR finder script and get out all exon span
	
0405:
	review about all the depth-related info of all missed exons
	prepare for a review presentation
	prepare for the LCR calling methods
	
	
0406:
	use rsync for transfer data instead of scp
	using UID for classify the mixed samples from the probing procedure
	
	also removing duplicates procedure should be performed
	
	DMD re-sequencing analysis:
		select out duplication samples and het deletions

	
0407:
	get out non-cordinate samples
	arrange for P17 and 2M sequencing for different samples
	
0410:
	keep in mind for the 2M samples, there are now 2 runs that are not currently analyzed
	choose out biased sample out and arrange for 2m AND p17 SEQUENCING
	
	prepare and get back to work
	objectif:
		1> get out diff sample and prepare for email
		2> get out quality control for 2M sequencing pipeline
		3> validate with SYSU-AF with biased sample
		4> optimize the LCR finder using LD algo
	
	python -m py_compile file.py
	
0411:
	prepare for the 2M QC panal
	including using depth plot and coverage
	
	note that recent P37 has lots of indel, should keep good QC control
	
0412:
	test the new 2M QC pipeline
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	