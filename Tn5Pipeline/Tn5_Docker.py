

##Tn5 analysis pipeline
import os
import time
import sys
import subprocess
from subprocess import call, Popen, PIPE, STDOUT
import csv
import io
import glob
from itertools import groupby
import re
import shutil

if __name__ == "__main__":
    master_folder_path = str(sys.argv[1])
    sampleID_fasta_csv = str(sys.argv[2])

if (master_folder_path[-1]=="/"):
	master_folder_path=master_folder_path[:-1]


##Need to add to reflow script
adapter_path = "NexteraPE-PE.fa"

# your trimmomatic parameters
lead_score=str(3)
trail_score=str(3)
min_len=str(50)
wind_size=str(4)
wind_qual=str(20)

def get_col(arr, col):
    return map(lambda x : x[col], arr)

def create_folder(fastq_folder_path):
	mapping = []
	with open(sampleID_fasta_csv, "r") as f:
		reader = csv.reader(f, delimiter=",")
		count = 0
		for row in reader:
			mapping.append(row)
	ID_fasta_list = mapping[1:]
	listed_sampleIDs = list(get_col(ID_fasta_list,0))
	print(list(set(listed_sampleIDs)))

	sampleID = []
	full_fq_path = []
	for fastq_file in glob.glob(master_folder_path+"/*"):
		samp_id = str(os.path.basename(fastq_file).split("_S")[0])
		print("sample id from fastq folder: " + (samp_id))
		if listed_sampleIDs.count(samp_id) > 0 :
			print("double cross check passed: " + str(os.path.basename(fastq_file.split("_S")[0])))
			full_fq_path.append(fastq_file)
			sampleID.append(fastq_file.split("_S")[0])
	sampleID_clean = list(set(sampleID)) # assumes all samples are PAIRED and therefore have two associated fastqs (removes duplicates)
	
	# create directory for every sample ID
	for directory in sampleID_clean:
		os.makedirs("/home/ubuntu/" + directory)

	# move fastq files to the corresponding sample ID directory
	for fq_path in full_fq_path:
		new_path = fq_path.split("_S")[0]+"/"+os.path.basename(fq_path)
		shutil.move(fq_path,new_path)
	
	# move fasta files to the corresponding sample ID directory
	for entry in ID_fasta_list:
		sample_id = entry[0]
		if os.path.exists("/home/ubuntu/" + master_folder_path + "/"+sample_id):
			fa_path = "/home/ubuntu/" + master_folder_path+"/"+sample_id+"/"+sample_id+".fasta"
			sequence = entry[1]
			f = open(fa_path,"w")
			f.write(">"+sample_id+"\n")
			f.write(sequence)
			f.close
		 
	
	return master_folder_path
	
sorted_master_path = create_folder(master_folder_path)


def run(folder_path):
	print("folder path: " + str(folder_path))
	print(os.path.join(folder_path, '*.fasta'))
	# extract fasta path and set as only key to dictionary (with no values either)
	for filename in glob.glob(os.path.join(folder_path+"/", '*.fasta')):
		fasta_path = os.path.basename(filename)
		print("fasta path is: " + str(fasta_path))
		d1 = {fasta_path:None}

	def natural_sort(l): 
	    convert = lambda text: int(text) if text.isdigit() else text.lower() 
	    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	    return sorted(l, key = alphanum_key)

	# extract fastq files and generate list of fastq file paths
	fastq_list = []
	for filename in glob.glob(os.path.join(folder_path, '*.fastq')):
		fastq_list.append(os.path.basename(filename))

	sorted_fastqs = zip(natural_sort(fastq_list)[0::2], natural_sort(fastq_list)[1::2])
	

	def dockerpulls():
		p = Popen(["docker", "pull", "fjukstad/trimmomatic"])
		p.communicate()
		p.wait()

		q = Popen(["docker", "pull", "fjukstad/bwa"])
		q.communicate()
		q.wait()

		r = Popen(["docker", "pull", "fjukstad/samtools"])
		r.communicate()
		r.wait()

		s = Popen(["docker", "pull", "biocontainers/bcftools"])
		s.communicate()
		s.wait()

		t = Popen(["docker", "pull", "fjukstad/picard/"])
		t.communicate()
		t.wait()


		return None

	def trim_filter_fqs(paired_fq_list):

		for entry in paired_fq_list:
			fwd_read = entry[0]
			rev_read = entry[1]
			details = fwd_read.split("_S")
			name = details[0]

			trimmed_R1_path = folder_path + "/"+ name + "_trimmed_R1.fastq"
			trimmed_R2_path = folder_path + "/"+ name + "_trimmed_R2.fastq"
			unpaired_R1_path = folder_path + "/"+ name + "_unpaired_R1.fastq"
			unpaired_R2_path = folder_path + "/"+ name + "_unpaired_R2.fastq"
			os.chdir(folder_path)
			outputs = [trimmed_R1_path,unpaired_R1_path,trimmed_R2_path,unpaired_R2_path]

			
			cut_folder_path = folder_path.split("/home/ubuntu/")
			cut_folder_path = cut_folder_path[1]

			cut_trimmed_R1_path = trimmed_R1_path.split("/home/ubuntu/")
			cut_trimmed_R1_path = cut_trimmed_R1_path[1]

			cut_trimmed_R2_path = trimmed_R2_path.split("/home/ubuntu/")
			cut_trimmed_R2_path = cut_trimmed_R2_path[1]

			cut_unpaired_R1_path = unpaired_R1_path.split("/home/ubuntu/")
			cut_unpaired_R1_path = cut_unpaired_R1_path[1]

			cut_unpaired_R2_path = unpaired_R2_path.split("/home/ubuntu/")
			cut_unpaired_R2_path = cut_unpaired_R2_path[1]

			p = Popen(["sudo","docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/trimmomatic","PE","-phred33",cut_folder_path+"/"+fwd_read,cut_folder_path+"/"+rev_read,
				cut_trimmed_R1_path,cut_unpaired_R1_path,cut_trimmed_R2_path,cut_unpaired_R2_path,"ILLUMINACLIP:"+adapter_path
				+":4:20:10","LEADING:"+lead_score,"TRAILING:"+trail_score,"SLIDINGWINDOW:"+wind_size+":"+wind_qual,
				"MINLEN:"+min_len],stdout=PIPE)
			p.communicate()
			p.wait()


		return None

	trim_filter_fqs(sorted_fastqs)

	# retrieve the trimmed, paired fastqs (throw out the unpaired fastqs)
	trimmed_fqs = []
	for filename in glob.glob(os.path.join(folder_path, '*_trimmed_*.fastq')):
		fq = os.path.basename(filename)
		trimmed_fqs.append(fq)

	sorted_trimmed_fqs = zip(natural_sort(trimmed_fqs)[0::2], natural_sort(trimmed_fqs)[1::2])

	d1[fasta_path] = sorted_trimmed_fqs


	def aln_to_ref(fa_fq_dict):

		for key in fa_fq_dict:
			ref_path = folder_path + "/" + key
			paired_fastqs = fa_fq_dict.get(key)
			for paired_read in paired_fastqs:
				name = (paired_read[0].split('_t'))
				# define fastq, sam, bam paths
				fq_fwd_path = folder_path + "/" + paired_read[0]
				fq_rev_path = folder_path + "/" + paired_read[1]
				output_sam = folder_path+ "/" + name[0] + ".sam"
				sorted_bam = folder_path+ "/" + name[0] + ".bam"

				
				#Use system to call BWA
				#index reference

				cut_ref_path = ref_path.split("/home/ubuntu/")
				cut_ref_path = cut_ref_path[1]

				cut_fq_fwd_path = fq_fwd_path.split("/home/ubuntu/")
				cut_fq_fwd_path = cut_fq_fwd_path[1]

				cut_fq_rev_path = fq_rev_path.split("/home/ubuntu/")
				cut_fq_rev_path = cut_fq_rev_path[1]

				cut_output_sam = output_sam.split("/home/ubuntu/")
				cut_output_sam = cut_output_sam[1]

				cut_sorted_bam = sorted_bam.split("/home/ubuntu/")
				cut_sorted_bam = cut_sorted_bam[1]

				process1 = call(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/bwa", "index",cut_ref_path], stdout = PIPE) # index reference

				# generate sam files
				with open(output_sam, "w") as f1:
					p = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/bwa","mem",cut_ref_path,cut_fq_fwd_path,cut_fq_rev_path],stdout=f1)
					p.communicate()
					p.wait()

				#Conversion of .sam to .bam
				with open(sorted_bam, "w") as f2:
					process2 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/samtools", "sort",cut_output_sam],stdout=f2)
					process2.communicate()
					process2.wait()
					process3 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/samtools", "index",cut_sorted_bam])
				
		return None
		

	aln_to_ref(d1)

	#######Dockerized to here
	

	bam_list = []

	for filename in glob.glob(os.path.join(folder_path, '*.bam')):
		bam_list.append(os.path.basename(filename))
	sorted_bam = natural_sort(bam_list)
 

	def samtools_flagstatoutput(bam_file_list):
	
		
		for entry in bam_file_list:
			flagstat_txt = open(folder_path+"/"+entry.split(".bam")[0]+"_flagstats.txt","w")

			cut_folder_path = folder_path.split("/home/ubuntu/")
			cut_folder_path = cut_folder_path[1]

			process3 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/samtools","flagstat",cut_folder_path+"/"+entry],stdout=PIPE)
			
			a = process3.stdout.read().decode("utf-8")
			flagstat_txt.write(entry.split(".bam")[0] + "\n")
			flagstat_txt.write("%s\n" % a)
		flagstat_txt.close()
		

		return None

	samtools_flagstatoutput(sorted_bam)


	def samtools_mpileup(bam_file_list):
		#SAMTOOLS mpileup
		ref_path = folder_path+"/"+fasta_path

		for entry in bam_file_list:
			#retrieve sample name
			details = entry.split(".bam")
			name = details[0]

			# generate sample-specific mpileup files
			bam_path = folder_path+"/"+entry
			mpileup_file = folder_path + "/" + name + ".mpileup"
			vcf_file = folder_path + "/" + name + ".vcf"
			stats_file = folder_path + "/" + name + "_" + "stats.txt"
			picard_wgs_metrics = folder_path + "/" + name + "_" + "picard_wgs.txt"
			picard_size_metrics = folder_path + "/" + name + "_" + "picard_size.txt"
			picard_size_hist = folder_path + "/" + name + "_" + "picard_size.pdf"
			parsed_picard_wgs_output = folder_path + "/" + name + "_" + "picard_wgs_parsed.csv"
			consensus = folder_path + "/" + name + "_" + "consensus.fq"
			
			cut_ref_path = ref_path.split("/home/ubuntu/")
			cut_ref_path = cut_ref_path[1]

			cut_bam_path = bam_path.split("/home/ubuntu/")
			cut_bam_path = cut_bam_path[1]

			# write mpileup file
			with open(mpileup_file, "w") as f2:
				process2 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","biocontainers/bcftools", "mpileup","-f",cut_ref_path, cut_bam_path],stdout=f2)
				process2.communicate()
				process2.wait()

			cut_mpileup_file = mpileup_file.split("/home/ubuntu/")
			cut_mpileup_file = cut_mpileup_file[1]

			# generate vcf files
			with open(vcf_file,"w") as f3:
				process3 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","biocontainers/bcftools","call","-c",cut_mpileup_file],stdout=f3)
				process3.communicate()
				process3.wait()
			
			cut_vcf_file = vcf_file.split("/home/ubuntu/")
			cut_vcf_file = cut_vcf_file[1]
			
			# generate stats
			with open(stats_file,"w") as f5:
				process5 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","biocontainers/bcftools","stats",cut_vcf_file],stdout=f5)
				process5.communicate()
				process5.wait()

			cut_picard_wgs_metrics = picard_wgs_metrics.split("/home/ubuntu/")
			cut_picard_wgs_metrics = cut_picard_wgs_metrics[1]


			# generate picard stats
			process6 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/picard","CollectWgsMetrics","COVERAGE_CAP=100000","SAMPLE_SIZE=5000","I="+cut_bam_path,"R="+cut_ref_path,"O="+picard_wgs_metrics],stdout=PIPE)
			process6.communicate()
			process6.wait()
			stats =[]
			with open(picard_wgs_metrics,'r') as f:
				reader=csv.reader(f,delimiter="\t")
				for line in reader:
					stats.append(line)
			stat_titles = stats[6]
			stat_values = stats[7]
			with open(parsed_picard_wgs_output,"w") as f2:
				writer=csv.writer(f2,delimiter=",")
				for i in range(len(stat_titles)):
					writer.writerow([stat_titles[i],stat_values[i]])

			cut_picard_size_hist = picard_size_hist.split("/home/ubuntu/")
			cut_picard_size_hist = cut_picard_size_hist[1]

			cut_picard_size_metrics = picard_size_metrics.split("/home/ubuntu/")
			cut_picard_size_metrics = cut_picard_size_metrics[1]

			# generate picard insert size stats
			process7 = Popen(["docker","run","-v","/home/ubuntu/:/DATA","-w","/DATA","fjukstad/picard","CollectInsertSizeMetrics","I="+cut_bam_path,"H="+cut_picard_size_hist,"O="+cut_picard_size_metrics],stdout=PIPE)
			process7.communicate()
			process7.wait()

			# generate consensus sequence
			with open(consensus, "w") as f9:
				p9 = Popen(["vcfutils.pl","vcf2fq", vcf_file],stdout=f9)
				p9.communicate()
				p9.wait()
			
		# mpileup_txt.close()
		return None


	samtools_mpileup(sorted_bam)

	return None


for sub_folder in glob.glob(master_folder_path+"/*"):
	print("sub_folder is :" + str(sub_folder))
	# IF SUB_FOLDER IS DIRECTORY THEN RUN!
	if sub_folder.count(".fastq") > 0:
		pass
	else:
		run("/home/ubuntu/" + sub_folder)

#no dockerization here
def compile_files(data_folder):

	for folder in glob.glob("/home/ubuntu/"+master_folder_path+"/*"):
		temp_size = []
		output_size_path = os.path.abspath(folder) +"/" + os.path.basename(folder)+"_"+ "size_distribution.csv"
		for file in glob.glob(folder+"/*"):
			if "_size.txt" in file:
				with open(file, "r") as f:
					reader=csv.reader(f,delimiter="\t")
					for line in reader:
						temp_size.append(line)
				size_stats=temp_size[11:-1]
				with open(output_size_path,"w") as f:
					writer = csv.writer(f,delimiter=",")
					for entry in size_stats:
						writer.writerow(entry)

	#compile and merge picard Wgs metric outputs
	wgs_metric_paths = []
	stat_txt_paths = []
	size_metric_paths = []
	for folder in glob.glob("/home/ubuntu/"+master_folder_path+"/*"):
		print("folder is: "+str(folder))
		for file in glob.glob(folder+"/*"):
			print("file is: "+str(file))
			if "parsed" in file:
				wgs_metric_paths.append(file)
			if "_stats.txt" in file:
				stat_txt_paths.append(file)
			if "_size.txt" in file:
				size_metric_paths.append(file)

	merged_stats_file = "/home/ubuntu/" + data_folder +"/merged_stats.csv"
	stats = [["Metric"]]
	# get picard wgs column names of interest
	n=0
	for item in wgs_metric_paths:
		while n<1:
			with open(item,"r") as f:
				reader = csv.reader(f,delimiter=",")
				for line in reader:
					stats.append([line[0]])
			n=1
	# get samtools stats column names of interest
	n=0
	temp_stat =[]
	for item in stat_txt_paths:
		while n<1:
			with open(item,"r") as f:
				reader = csv.reader(f,delimiter="\t")
				for line in reader:
					temp_stat.append(line)
			n=1

	# get picard size column names
	n=0
	temp_stat3 =[]
	for item in size_metric_paths:
		while n<1:
			with open(item,"r") as f:
				reader = csv.reader(f,delimiter="\t")
				for line in reader:
					temp_stat3.append(line)
			n=1

	# add samtools stats columns to master stat columns
	stats2 = temp_stat[25:31]
	for x in stats2:
		stats.append([x[2]])
	

	# add picard size columns to master stat columns
	stats3_titles = temp_stat3[6][:7]
	for item in stats3_titles:
		stats.append([item])

	# gather picard wgs stats and append to master stats list
	for item in wgs_metric_paths:
		n=0
		with open(item,"r") as f:
			reader = csv.reader(f,delimiter=",")
			name = os.path.basename(item).split("_p")[0]
			stats[n].append(name)
			for line in reader:
				n+=1
				stats[n].append(line[1])
	
	# gather samtools stats and append to master stats list
	for item in stat_txt_paths:
		temp_stat2=[]
		n=0
		with open(item,"r") as f:
			reader = csv.reader(f,delimiter="\t")
			for line in reader:
				temp_stat2.append(line)
			x = temp_stat2[25:31]
			n=29
			for item in x:
				stats[n].append(item[3])
				n+=1
	
	# gather picard size stats and append to master stats list
	for item in size_metric_paths:
		temp_stat2=[]
		n=0
		with open(item,"r") as f:
			reader = csv.reader(f,delimiter="\t")
			for line in reader:
				temp_stat2.append(line)
			x = temp_stat2[7][:7]
			n=35
			for item in x:
				stats[n].append(item)
				n+=1
	
	print("filled stat list: " + str(stats))
	with open(merged_stats_file,"w") as f:
		writer = csv.writer(f,delimiter=",")
		for item in stats:
			writer.writerow(item)
				
	with open(merged_stats_file,"w") as f:
		writer = csv.writer(f,delimiter=",")
		for item in stats:
			writer.writerow(item)

	output_path = "/home/ubuntu/"+data_folder+"/output/"
	os.makedirs(output_path)
	move_keys = [".bam", "consensus", "merged_stats","size.pdf"]

	for folder in glob.glob("/home/ubuntu/"+master_folder_path+"/*"):
		if any(x in os.path.basename(folder) for x in move_keys):
			old = os.path.abspath(folder)
			new = output_path+os.path.basename(folder)
			shutil.move(old,new)
		for file in glob.glob(folder+"/*"):
			if any(x in os.path.basename(file) for x in move_keys):
				old = os.path.abspath(file)
				new = output_path+os.path.basename(file)
				shutil.move(old,new)

	return None

	#rm non-called fastq in current (but not subdirectories)
	#rename output - summary stats


compile_files(master_folder_path)
