import re

def fasta_to_dict(filename):
	"""
	>chr1
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	"""
	outdict = {}
	chromlist = []
	f = open(filename,'r')
	for str_x in f:
		str_x = str_x.strip("\n")
		if str_x[0] == ">":
			chrom = str_x[1:]
			chromlist.append(chrom)
			outdict[chrom] = []
		else:
			outdict[chrom].append(str_x)
	chromlist.append(chrom)
	print("step1 is over")
	for chrom in chromlist:
		sublist = outdict[chrom]
		substr = "".join(sublist)
		outdict[chrom] = substr
	print("step2 is over")
	return outdict

def single_base_method_bed_to_dict(single_base_method_bed_filename):
	"""
	chr1    564544  564545  chr1_564544_564545_AAACA_A_f
	chr1    564559  564560  chr1_564559_564560_AAACA_A_f
	chr1    565164  565165  chr1_565164_565165_AAACC_A_f
	"""
	f = open(single_base_method_bed_filename,'r')
	miCLIP_dict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		start = int(list_x[1])
		end = int(list_x[2])
		sitename = list_x[3]
		if sitename.split("_")[-1] == "f":
			strand = "+"
		elif sitename.split("_")[-1] == "r":
			strand = "-"
		else:
			raise Exception("miCLIP file strand error")
		try:
			miCLIP_dict[(chrom,strand)].append([chrom,start,end,sitename,strand])
		except:
			miCLIP_dict[(chrom,strand)] = [[chrom,start,end,sitename,strand]]
	return miCLIP_dict

def single_base_method_overlap(chrom,start,end,strand,miCLIP_dict):
	try:
		miCLIP_list = miCLIP_dict[(chrom,strand)]
	except:
		return "False"
	for miCLIP in miCLIP_list:
		miStart =  miCLIP[1]
		miEnd = miCLIP[2]
		miName = miCLIP[3]
		if miStart < end and miStart >= start:
			return "True"
	return "False"
	

def peak_format_fileread(filename,fasta_dict,bed_base_type,miCLIP_dict,DARTseq_dict,result_file):
	"""
	chrom   start   end     peakname						strand  
	chr12   6980098 6980198 ENSG00000111671_exon_1_1        -  
	chr12   6980148 6980248 ENSG00000111671_exon_1_2        -      
	"""
	f = open(filename,'r')
	d = open(result_file,'a')
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if str_x[0]=="#" or list_x[0]=="chrom":
			d.write("\t".join(["chrom","start","end","peakname","strand","outsequence","match_GGACT","match_RRACH","match_DRACH","match_miCLIP","match_DARTseq"])+"\n")
			continue
		chrom = list_x[0]
		start = int(list_x[1])
		end = int(list_x[2])
		peakname = list_x[3]
		strand = list_x[4]
		if bed_base_type == "0base":
			sequence = fasta_dict[chrom][start:end]
		if bed_base_type == "1base":
			sequence = fasta_dict[chrom][start-1:end]
		outsequence = sequence_reverse(sequence=sequence,strand=strand)
		match_GGACT = motif_parten_match(sequence=outsequence,motif_parten="GGACT")
		match_RRACH = motif_parten_match(sequence=outsequence,motif_parten="[GA][GA]AC[ATC]")
		match_DRACH = motif_parten_match(sequence=outsequence,motif_parten="[GAT][GA]AC[ATC]")
		if bed_base_type == "0base":
			match_miCLIP = single_base_method_overlap(chrom=chrom,start=start,end=end,strand=strand,miCLIP_dict=miCLIP_dict)
			match_DARTseq = single_base_method_overlap(chrom=chrom,start=start,end=end,strand=strand,miCLIP_dict=DARTseq_dict)

		if bed_base_type == "1base":
			match_miCLIP = single_base_method_overlap(chrom=chrom,start=start-1,end=end,strand=strand,miCLIP_dict=miCLIP_dict)
			match_DARTseq = single_base_method_overlap(chrom=chrom,start=start-1,end=end,strand=strand,miCLIP_dict=DARTseq_dict)

		str_write = "\t".join([chrom,str(start),str(end),peakname,strand,outsequence,match_GGACT,match_RRACH,match_DRACH,match_miCLIP,match_DARTseq]) + "\n"
		d.write(str_write)
	d.close()
	f.close()

def sequence_reverse(sequence,strand):
	reverse_dict = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
	sequence  = sequence.upper()
	if strand == "+":
		outsequence = sequence
	if strand == "-":
		outlist = []
		for base in sequence:
			reverse_base = reverse_dict[base]
			outlist.append(reverse_base)
		outlist.reverse()
		outsequence = "".join(outlist)
	return outsequence

def motif_parten_match(sequence,motif_parten):
	motif_list = re.findall(motif_parten,sequence)
	if motif_list == []:
		return "False"
	else:
		return "True"


def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='m6A peak quality control by motif and miCLIP dataset overlap')
	parser.add_argument('--m6A_peak_bed_format_filename', required=True, help="m6A peak bed format filename")
	parser.add_argument('--fasta_USCS_download_filename', required=True, help="fasta USCS download filename")
	parser.add_argument('--miCLIP_dataset_filename', required=True, help="miCLIP dataset filename")
	parser.add_argument('--DARTseq_dataset_filename', required=True, help="miCLIP dataset filename")
	parser.add_argument('--base_type', required=True, help="base type")
	parser.add_argument('--result_filename', required=True, help="result filename")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	m6A_peak_bed_format_filename = args.m6A_peak_bed_format_filename
	fasta_USCS_download_filename = args.fasta_USCS_download_filename
	miCLIP_dataset_filename = args.miCLIP_dataset_filename
	DARTseq_dataset_filename = args.DARTseq_dataset_filename
	base_type = args.base_type
	result_filename = args.result_filename
	fasta_dict = fasta_to_dict(filename=fasta_USCS_download_filename)
	miCLIP_dict = single_base_method_bed_to_dict(single_base_method_bed_filename=miCLIP_dataset_filename)
	DARTseq_dict = single_base_method_bed_to_dict(single_base_method_bed_filename=DARTseq_dataset_filename)
	peak_format_fileread(filename=m6A_peak_bed_format_filename,
						 fasta_dict=fasta_dict,
						 bed_base_type=base_type,
						 miCLIP_dict=miCLIP_dict,
						 DARTseq_dict=DARTseq_dict,
						 result_file=result_filename)


if __name__=="__main__":
	main()



