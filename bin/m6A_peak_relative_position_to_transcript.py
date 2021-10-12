
def transcript_to_dict_bin_to_trans(filename,binSize,empty_dict):
	'''
	76      NM_001011874.1  chr1    -       3214481 3671498 3216021 3671348 3       3214481,3421701,3670551,        3216968,3421901,3671498,      0       Xkr4    cmpl    cmpl
	76      XM_011238395.2  chr1    -       3322815 3672278 3323746 3671348 3       3322815,3421701,3670551,        3323760,3421901,3672278,      0       Xkr4    cmpl    cmpl
	'''
	f = open(filename,'r')
	#outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if str_x[0]=="#":
			continue
		transid = list_x[1]
		chrom = list_x[2]
		strand = list_x[3]
		gene_start = int(list_x[4])
		gene_end = int(list_x[5])
		cds_start = int(list_x[6])
		cds_end = int(list_x[7])
		exon_number  = list_x[8]
		exon_start_list = [int(x) for x in list_x[9].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[10].split(",")[:-1]]
		geneid = list_x[11]
		min_binNum = gene_start//binSize
		max_binNum = gene_end//binSize
		
		for binNum in range(min_binNum,max_binNum+1):
			try:
				empty_dict[(chrom,binNum)].append([gene_start,gene_end,cds_start,cds_end,exon_start_list,exon_end_list,strand,geneid])
			except:
				pass
	return empty_dict

def transcript_to_dict_gene_to_trans(filename):
	'''
	76      NM_001011874.1  chr1    -       3214481 3671498 3216021 3671348 3       3214481,3421701,3670551,        3216968,3421901,3671498,	0	Xkr4    cmpl    cmpl
	76      XM_011238395.2  chr1    -       3322815 3672278 3323746 3671348 3       3322815,3421701,3670551,        3323760,3421901,3672278,	0	Xkr4    cmpl    cmpl
	'''
	f = open(filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		if str_x[0]=="#":
			continue
		transid = list_x[1]
		chrom = list_x[2]
		strand = list_x[3]
		gene_start = int(list_x[4])
		gene_end = int(list_x[5])
		cds_start = int(list_x[6])
		cds_end = int(list_x[7])
		exon_number  = list_x[8]
		exon_start_list = [int(x) for x in list_x[9].split(",")[:-1]]
		exon_end_list = [int(x) for x in list_x[10].split(",")[:-1]]
		geneid = list_x[11]
		outdict[geneid] = [gene_start,gene_end,cds_start,cds_end,exon_start_list,exon_end_list,strand]
	return outdict

def empty_dict_make(genome_size_filename,binSize):
	"""
	chr1    249250621
	chr10   135534747
	chr11   135006516
	chr12   133851895
	"""
	f = open(genome_size_filename,'r')
	outdict = {}
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		genome_size = int(list_x[1])
		max_binNum = genome_size//binSize
		for binNum in range(max_binNum+1):
			outdict[(chrom,binNum)] = []
	return outdict

def bed_ovarlap_transcript(transDict,chrom,start,end,strand,binSize):
	start_binNum = start//binSize
	end_binNum = end//binSize
	translist = []
	for binNum in range(start_binNum,end_binNum+1):
		translist += transDict[(chrom,binNum)]
	geneid_list = []
	for sublist in translist:
		gene_start,gene_end,cds_start,cds_end,exon_start_list,exon_end_list,gene_strand,geneid = sublist
		if gene_strand == strand:
			overlap_or_not = overlap_transcript(chrom=chrom,start=start,end=end,exon_start_list=exon_start_list,exon_end_list=exon_end_list)
			if overlap_or_not == True:
				geneid_list.append(geneid)
			else:
				pass
		else:
			pass
	return geneid_list

def overlap_transcript(chrom,start,end,exon_start_list,exon_end_list):
	chromStart = start + (end - start)//2
	chromEnd = chromStart + 1
	for i in range(len(exon_start_list)):
		exon_start = exon_start_list[i]
		exon_end = exon_end_list[i]
		if chromStart < exon_end and chromStart >= exon_start:
			return True
		else:
			pass
	return False
	

def precentage_for_transcript_part(trans_geneStart,trans_geneEnd,trans_cdsStart,trans_cdsEnd,trans_chromStart,trans_chromEnd,transcript_strand):
	UTR3_start = trans_geneStart
	UTR3_end = trans_cdsStart
	cds_start = trans_cdsStart
	cds_end = trans_cdsEnd
	UTR5_start = trans_cdsEnd
	UTR5_end = trans_geneEnd
	if (UTR3_end - UTR3_start)<=0 or (cds_end - cds_start)<=0 or (UTR5_end - UTR5_start)<=0:
		return -1
	overlap_type1,overlap_start1,overlap_end1,overlap_len1 = bed_overlap_bed(trans_chromStart,trans_chromEnd,UTR3_start,UTR3_end)
	overlap_type2,overlap_start2,overlap_end2,overlap_len2 = bed_overlap_bed(trans_chromStart,trans_chromEnd,cds_start,cds_end) 
	overlap_type3,overlap_start3,overlap_end3,overlap_len3 = bed_overlap_bed(trans_chromStart,trans_chromEnd,UTR5_start,UTR5_end)
	if overlap_len1 >=1:
		p = int(30*(trans_chromStart - UTR3_start)/(UTR3_end - UTR3_start))
	elif overlap_len2 >=1:
		p = int(40*(trans_chromStart - cds_start)/(cds_end - cds_start)) + 30
	elif overlap_len3 >=1:
		p = int(30*(trans_chromStart - UTR5_start)/(UTR5_end - UTR5_start)) + 70
	else:
		return -1
		#raise Exception("precentage_for_transcript_part error")
	if transcript_strand == "+":
		return p
	elif transcript_strand == "-":
		return 99 - p
	else:
		raise Exception("precentage_for_transcript_part error")

def genomePosition_to_transcriptPosition(genome_position,exon_list):
	transcript_position = 0
	for exon in exon_list:
		exon_start = exon[0]
		exon_end = exon[1]
		if genome_position >= exon_end and genome_position > exon_start:
			transcript_position += (exon_end - exon_start)
		elif genome_position < exon_end and genome_position >= exon_start:
			transcript_position += (genome_position - exon_start)
		elif genome_position < exon_end and genome_position < exon_start:
			transcript_position += 0
		else:
			raise Exception("genomePosition_to_transcriptPosition error")
	return transcript_position

def exon_length_calculate(exon_list):
	length = 0
	for i in range(len(exon_list)):
		exon_start_i = exon_list[i][0]
		exon_end_i = exon_list[i][1]
		length += (exon_end_i - exon_start_i)
	return length


def transcriptPosition_to_genomePosition(transcript_position,exon_list):
	genome_position = exon_list[0][0]
	tmp_transcript_position = transcript_position
	for i in range(len(exon_list)):
		exon_start = exon_list[i][0]
		exon_end = exon_list[i][1]
		tmp_genome_position = genome_position + tmp_transcript_position
		if tmp_genome_position > exon_end and tmp_genome_position > exon_start:
			genome_position = exon_list[i+1][0]
			tmp_transcript_position = tmp_transcript_position - (exon_end - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position > exon_start:
			genome_position = tmp_genome_position
			tmp_transcript_position = tmp_transcript_position - (tmp_genome_position - exon_start)
		elif tmp_genome_position == exon_end and tmp_genome_position > exon_start:
			if i == (len(exon_list)-1):
				genome_position = tmp_genome_position
				return genome_position
			else:
				genome_position = exon_list[i+1][0]
				tmp_transcript_position = tmp_transcript_position - (exon_end - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position == exon_start:
			genome_position = tmp_genome_position
			tmp_transcript_position = tmp_transcript_position - (tmp_genome_position - exon_start)
		elif tmp_genome_position < exon_end and tmp_genome_position < exon_start:
			pass
		else:
			raise Exception("transcriptPosition_to_genomePosition error")
	return genome_position



def bed_overlap_bed(start1,end1,start2,end2):
	"""
	input:1000,1200,1100,1300
	output:5,1100,1200,100
	use:overlap_type,overlap_start,overlap_end,overlap_len = bed_overlap_bed(start1,end1,start2,end2)
	"""
	########################################################################          1
	if start1 > start2 and start1 > end2 and end1 > start2 and end1 > end2:
		overlap_type = 1
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	# type1 equal
	elif start1 > start2 and start1 == end2 and end1 > start2 and end1 > end2:
		overlap_type = 1
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	########################################################################          2
	elif start1 > start2 and start1 < end2 and end1 > start2 and end1 > end2:
		overlap_type = 2
		overlap_start = start1
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	########################################################################          3
	elif start1 < start2 and start1 < end2 and end1 > start2 and end1 > end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	# type3 equal
	elif start1 == start2 and start1 < end2 and end1 > start2 and end1 > end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	elif start1 < start2 and start1 < end2 and end1 > start2 and end1 == end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	elif start1 == start2 and start1 < end2 and end1 > start2 and end1 == end2:
		overlap_type = 3
		overlap_start = start2
		overlap_end = end2
		overlap_len = overlap_end - overlap_start
	#########################################################################         4
	elif start1 > start2 and start1 < end2 and end1 > start2 and end1 < end2:
		overlap_type = 4
		overlap_start = start1
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	# type4 equal
	elif start1 == start2 and start1 < end2 and end1 > start2 and end1 < end2:
		overlap_type = 4
		overlap_start = start1
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	elif start1 > start2 and start1 < end2 and end1 > start2 and end1 == end2:
		overlap_type = 4
		overlap_start = start1
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	#########################################################################          5
	elif start1 < start2 and start1 < end2 and end1 > start2 and end1 < end2:
		overlap_type = 5
		overlap_start = start2
		overlap_end = end1
		overlap_len = overlap_end - overlap_start
	#########################################################################          6
	elif start1 < start2 and start1 < end2 and end1 < start2 and end1 < end2:
		overlap_type = 6
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	# type6 equal
	elif start1 < start2 and start1 < end2 and end1 == start2 and end1 < end2:
		overlap_type = 6
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	########################################################################           7
	elif start1 == end1 or start2 == end2:
		overlap_type = 7
		overlap_start = -1
		overlap_end = -1
		overlap_len = overlap_end - overlap_start
	else:
		raise Exception("bad overlap")
	return overlap_type,overlap_start,overlap_end,overlap_len


def peak_file_read(filename,max_length_transcript_dict,transdict_bin2trans,binSize,already_know_geneid,output_file):
	'''
	chr10   3957945 3958045 Plekhg1_exon_21_29      +
	chr10   3963500 3963600 Plekhg1_exon_22_3       +
	chr10   3963900 3964000 Plekhg1_exon_22_11      +
	'''
	f = open(filename,'r')
	d = open(output_file,'a')
	total_trans_p = []
	for str_x in f:
		str_x = str_x.strip("\n")
		list_x = str_x.split("\t")
		chrom = list_x[0]
		start = int(list_x[1])
		end = int(list_x[2])
		peakname = list_x[3]
		#geneid = peakname.split("_")[0]
		strand = list_x[4]
		if already_know_geneid == "True":
			geneid = peakname.split("_")[0]
			overlap_geneid_list = [geneid]
		else:
			overlap_geneid_list = bed_ovarlap_transcript(transDict=transdict_bin2trans,chrom=chrom,start=start,end=end,strand=strand,binSize=binSize)
		trans_p_list = []
		for geneid in overlap_geneid_list:
			gene_start,gene_end,cds_start,cds_end,exon_start_list,exon_end_list,transcript_strand = max_length_transcript_dict[geneid]
			exon_list = [[exon_start_list[i],exon_end_list[i]] for i in range(len(exon_start_list))]
			exon_length = exon_length_calculate(exon_list)
			chromStart = start + (end - start)//2
			chromEnd = chromStart + 1 
			trans_chromStart = genomePosition_to_transcriptPosition(chromStart,exon_list)
			trans_chromEnd = genomePosition_to_transcriptPosition(chromEnd,exon_list)
			trans_cdsStart = genomePosition_to_transcriptPosition(int(cds_start),exon_list)
			trans_cdsEnd = genomePosition_to_transcriptPosition(int(cds_end),exon_list)
			trans_geneStart = genomePosition_to_transcriptPosition(int(gene_start),exon_list)
			trans_geneEnd = genomePosition_to_transcriptPosition(int(gene_end),exon_list)
			trans_p = precentage_for_transcript_part(trans_geneStart,trans_geneEnd,trans_cdsStart,trans_cdsEnd,trans_chromStart,trans_chromEnd,transcript_strand)
			if trans_p == -1:
				pass
			else:
				trans_p_list.append(trans_p)
		trans_p_list = [x for x in trans_p_list if x!=-1]
		if trans_p_list==[]:
			trans_p = -1
			str_write = str_x+"\t"+str(trans_p) + "\n"
		else:
			trans_p = int(sum(trans_p_list)/len(trans_p_list))
			str_write = str_x+"\t"+str(trans_p) + "\n"
		total_trans_p.append(trans_p)
		#print(str_write)
		d.write(str_write)
	f.close()
	d.close()
	return total_trans_p

def precentage_result_write(total_trans_p_list,output_file):
	d = open(output_file+".precentage_for_plot.txt",'a')
	d.write("\t".join(["relative_location","counts"])+"\n")
	count_dict = {}
	for x in range(-1,100):
		count_dict[x] = 0
	for trans_p in total_trans_p_list:
		count_dict[trans_p] += 1
	for x in range(0,100):
		count = count_dict[x]
		str_write = "\t".join([str(x+1),str(count)]) + "\n"
		d.write(str_write)
	d.close()

def make_args():
	import argparse
	parser = argparse.ArgumentParser(description='peaks in transcript relative position')
	parser.add_argument('--peak_format_file', required=True, help="peak format file")
	parser.add_argument('--max_length_transcript_file', required=True, help="max length transcript file")
	parser.add_argument('--genome_size_file', required=True, help="genome size file")
	parser.add_argument('--binSize', required=True, help="binSize:100000")
	parser.add_argument('--already_know_geneid', required=True, help="already_know_geneid:True or False")
	parser.add_argument('--output_file', required=True, help="result output")
	args = parser.parse_args()
	return args

def main():
	args = make_args()
	peak_format_file = args.peak_format_file
	max_length_transcript_file = args.max_length_transcript_file
	genome_size_file = args.genome_size_file
	binSize = int(args.binSize)
	already_know_geneid = args.already_know_geneid
	output_file = args.output_file

	empty_dict = empty_dict_make(genome_size_filename=genome_size_file,binSize=binSize)
	transdict_bin2trans = transcript_to_dict_bin_to_trans(filename=max_length_transcript_file,binSize=binSize,empty_dict=empty_dict)
	transdict_geneid2trans = transcript_to_dict_gene_to_trans(filename=max_length_transcript_file)
	total_trans_p_list = peak_file_read(filename=peak_format_file,max_length_transcript_dict=transdict_geneid2trans,transdict_bin2trans=transdict_bin2trans,binSize=binSize,already_know_geneid=already_know_geneid,output_file=output_file)
	precentage_result_write(total_trans_p_list=total_trans_p_list,output_file=output_file)

if __name__=="__main__":
	main()



