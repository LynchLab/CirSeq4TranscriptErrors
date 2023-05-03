from scipy.stats import binom
import sys,re

#!!!!!
#1.call variants only if this variation is supported by one RNA fragment. 
#2.Skip indels in this script.
#3.This version allows the multi-category substitutions at one sites.
#4.Count substitutions based on reads instead of genomic sites.
#5.add filters according to the comparison between the expected binomial frequency and the observed frequency of one variant.
#6. add criteria on the minimum coverage for one locus

L=['+','-']
indel_count=0
total_coverage=0
cutoff=float(sys.argv[2])  #starts from 0.5
min_cov=int(sys.argv[3])

pos_cov,pos_sub,pos_ref={},{},{}	
base_cov={}	#coverages for four bases

##1st round of analysis
file_P=open(sys.argv[1],'r')
for each_line in file_P:
	line=each_line.split()
	if int(line[3])>0:
		comma,un,nun=0,0,0
		count_dic={'A':0,'T':0,'C':0,'G':0,'a':0,'t':0,'c':0,'g':0}
		sum_dic={}
		bases=line[4]
		total=int(line[3])
		ref=line[2]
		if any(x in bases for x in L):
			indel_count+=1
		else:	#Only visit substitutions
			key_potential=line[0]+"&&&"+str(line[1])
			for i in range(len(bases)): 
				if bases[i] == ',':
					comma+=1
				elif bases[i] in ('A','T','C','G','a','t','c','g'):
					count_dic[bases[i]]+=1
				elif bases[i] in ('N','n'):
					un+=1
				elif bases[i] == '^': # ^ marks the start of a read segment and is followed by quality base; $ marks the end of a read segment and no quality base follows. 
					if bases[i+1] in ('A','T','C','G','a','t','c','g'):
						count_dic[bases[i+1]]-=1
					if bases[i+1] in ('N','n'):
						nun+=1
			Un=un-nun
			total-=Un  # exclude Ns
			if total >= min_cov:  #? crucial to be >= 
				Mut=sum(iter(count_dic.values()))
				sum_dic['A']=count_dic['A']+count_dic['a'] 
				sum_dic['T']=count_dic['T']+count_dic['t']
				sum_dic['C']=count_dic['C']+count_dic['c']
				sum_dic['G']=count_dic['G']+count_dic['g']
				num_types=0
				mut={}
				pos_cov[key_potential]=total	#keep track on the effective coverage at each position
				pos_ref[key_potential]=ref	#keep track on the ref base ata each position
				for key in sum_dic:
					if sum_dic[key]>0: #Bases are denoted as "," or "." if read calls are consistant with reference. 
						num_types+=1
						mut[key]=sum_dic[key]
				pos_sub[key_potential]=[]	#if no errors, empty
				if 0< Mut/total< cutoff:	# ? <= or <
					for key in mut:
						sub=line[2]+':'+key
						for i in range(int(mut[key])):
							pos_sub[key_potential].append(sub)
				if ref in base_cov:	#keep track on coverages of different bases
					base_cov[ref]+=total
				else:
					base_cov[ref]=total
			 
##get estimate of miu from 1st round analysis
i,j=0,0
for key in pos_cov:
	i+=pos_cov[key]
	j+=len(pos_sub[key])
miu_1=j/i
#print (n,m,miu_1)

##iterations to filter out outliers 
miu_old=miu_1
miu_new=0
miu=[miu_1,0]
sub_count=[j]
B=len(pos_cov)
n=1

while miu_new < miu_old: # ? !=
	cov,count=0,0
	del_key=[]
	for key in pos_cov:
		x=len(pos_sub[key]) # the # of errors at one site
		y=pos_cov[key]	# the coverage of one site
		prob=1-binom.cdf(x-1, y, miu_old)	# prob of observe >= x 
		P_corrected=0.05/B	#bonferroni corrected P, use the number of evaluated loci instead of the number of loci with observed errors
		if prob < P_corrected:
			del_key.append(key)
		else:
			cov+=pos_cov[key]
			count+=len(pos_sub[key])
	if count == 0:	#for extreme case
		print ("Warning: no errors passed the filter")
		break
	for i in del_key:
		ref_del=pos_ref[i]
		base_cov[ref_del]-=pos_cov[i]
		del pos_cov[i]
		del pos_sub[i]
	miu_tmp=count/cov
	miu[n]=miu_tmp
	sub_count.append(count)
	miu_new=miu[n]
	miu_old=miu[n-1]
	B=len(pos_cov)
	n+=1
	miu.append(0)	#for miu_tmp in the next round 
	if n == 15:
		break

#final output
outfile1=open(sys.argv[1]+"_"+"analysis_cov50",'w')	# error info
outfile2=open(sys.argv[1]+"_"+"analysis_summary_cov50",'w')	#info on total coverages and iterations

#output1
outfile1.write('Chr'+'\t'+'Pos'+'\t'+'Sub'+'\t'+'Effective coverages'+'\t'+'Number of errors'+'\t'+'Frequency of errors'+'\n')
for key in pos_cov:
	info=key.split('&&&')
	if len(pos_sub[key]) > 0:
		sub_type={}	# there may be multiple types of sub at one position
		for k in pos_sub[key]:	 
			if k in sub_type:
				sub_type[k]+=1
			else:
				sub_type[k]=1
		for l in sub_type:
			for m in range(sub_type[l]): 
				outfile1.write(info[0]+'\t'+str(info[1])+'\t'+l+'\t'+str(pos_cov[key])+'\t'+str(sub_type[l])+'\t'+str(sub_type[l]/pos_cov[key])+'\n')

#output2
outfile2.write(str(miu)+'\n')
outfile2.write(str(sub_count)+'\n')
outfile2.write("Effective coverages:"+str(cov)+'\n')
for key in base_cov:
	outfile2.write(key+':'+str(base_cov[key])+'\n')
