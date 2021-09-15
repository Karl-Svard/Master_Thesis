import re
import subprocess

n=100
t_nom=1
t_den=100
out_name='n_'+str(n)+'_Tnom_'+str(t_nom)+'_Tden_'+str(t_den)+'.txt'
the_prec=1000
OUT_DICT={}
comm_list=['./distr_tav84']
comm_list+=[str(n)]
comm_list+=[str(t_nom)]
comm_list+=[str(t_den)]
comm_list+=[str(the_prec)]
#print comm_list
#raw_input()
a_proc=subprocess.Popen(comm_list, shell=False, stdin=subprocess.PIPE,stdout=subprocess.PIPE)
the_lines=re.split('\n',a_proc.communicate()[0].strip())
sum_up=0
outf=open('DIR_tav84_res/'+out_name,'w')
for x in the_lines:
	data=re.split('\s',x)
	sum_up+=float(data[1])
	print data
	out_str=data[0]+' '+data[1]
	outf.write(out_str+'\n')
print 'tot prob',sum_up
