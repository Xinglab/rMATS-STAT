#This script calculates the normalized posterior probability

import re,os,sys,rpy,math

rpy.r.set_seed(854920324);
rpy.r.source("Simu_Outlier.R");

#variation
sigma=math.sqrt(float(sys.argv[2]));
if len(sys.argv)>=4:
	min=float(sys.argv[3]);
	max=float(sys.argv[4]);
sigma2=math.sqrt(float(sys.argv[5]));

#ESRP_OE_RES.txt
ifile=open(sys.argv[1]);
ifile.readline();
ilines=ifile.readlines();
for i in ilines:
	K=5;
	element=re.findall('[^\t\n]+',i);
	normal_data=rpy.r.mysimu_data_unifplusnormal(float(element[0])+float(element[1]),float(element[2])+float(element[3]),sigma,min,max,sigma2)
	print_text='';#print('2');print(normal_data);
	for j in range(int(len(normal_data)/K)):
		this_print_text='';
		for k in range(K):
			this_print_text+=str(int(normal_data[k+K*j]))+',';
		#print(this_print_text);
		print_text+=this_print_text[:-1]+'\t';
	print(print_text[:-1]+'\t1\t1');
