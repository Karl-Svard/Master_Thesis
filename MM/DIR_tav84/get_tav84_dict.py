import re
import subprocess


def get_tav_dict(n,t_nom,t_den,the_prec):
	B_DICT={'tot prob':0.0,'vals':{0:0.0}}
	comm_list=['./DIR_tav84/distr_tav84']
	#comm_list=['/home/pers/DIR_divergence_estimate/DIR_MM/DIR_do_calculations/DIR_tav84/distr_tav84']
        #comm_list=['/home/pers/DIR_jobb/DIR_divergence_estimate/MM/DIR_do_calculations/DIR_tav84/distr_tav84']


	comm_list+=[str(n)]
	comm_list+=[str(t_nom)]
	comm_list+=[str(t_den)]
	comm_list+=[str(the_prec)]
	#print comm_list
	#raw_input()
	a_proc=subprocess.Popen(comm_list, shell=False, stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	the_lines=re.split('\n',a_proc.communicate()[0].strip().decode())
	sum_up=0
	for x in the_lines:
		data=re.split('\s',x)
		B_DICT['tot prob']+=float(data[1])
		B_DICT['vals'].update({int(data[0]):float(data[1])})
	return B_DICT

	
if __name__=="__main__":
	n=100
	t_nom=1
	t_den=100
	out_name='n_'+str(n)+'_Tnom_'+str(t_nom)+'_Tden_'+str(t_den)+'.txt'
	the_prec=1000
	d=get_tav_dict(n,t_nom,t_den,the_prec)
	for x in sorted(d['vals'].keys()):
		print (x,d['vals'][x])
	print ('tot prob',d['tot prob']) 
