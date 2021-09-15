import sys
import scipy.optimize
import scipy.integrate as integrate
import numpy as np

#sys.path.append("/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project/DIR_Per/method_scripts/mm/DIR_tav84/")
sys.path.append("/crex/proj/snic2019-8-14/private/Karl/MM/DIR_tav84")
#sys.path.append("/home/kalle/Documents/masters_project/uppmax_version/MM/DIR_tav84")
from get_tav84_dict import *


def get_prec_den():
	the_prec=1000
	tau_den=int(1e8)
	return [the_prec,tau_den]


def get_tau_num(a_tau,a_den):
	return int(round(a_tau*a_den))




def rising_fact_ratio(a,b,i):
	if i==0:
		return 1.0
	else:
		return 1.0*((1.0*a+i-1)/(1.0*b+i-1))*rising_fact_ratio(a,b,i-1)

def falling_fact_ratio(a,b,i):
	if i==0:
		return 1.0
	else:
		return 1.0*((1.0*a-(i-1))/(1.0*b-(i-1)))*falling_fact_ratio(a,b,i-1)


def h(m,l,i,j):
	if m==0 or l==0:
		if m==l:
			if i==0 and i==j:
				return 1.0
			return 0.0
		if m==0:
			if i==m and j>0:
				return 1.0
			return 0.0
		if j==l and i>0:
			return 1.0
		return 0.0
	if i==0 or j==0:
		return 0.0
	fact1=1.0*i*(m+l)/(1.0*m*l)
	fact2=rising_fact_ratio(j,1,i)
	fact3=falling_fact_ratio(m,m+l,i)
	fact4=falling_fact_ratio(l,m+l-i,j)
	return fact1*fact2*fact3*fact4

def h_star(m,l,i,j):
	if m==0 or l==0:
		if m==l:
			if i==0 and i==j:
				return 1.0
			return 0.0
		if m==0: #l is then >0
			if i==m and j>0:
				return 1.0
			return 0.0
		if j==l and i>0: #l is then ==0
			return 1.0
		return 0.0
	if i==0 or j==0:
		return 0.0
	fact1=1.0*i*(j+1)*(m+l)/(1.0*j*m*l)
	fact2=rising_fact_ratio(j,1,i)
	fact3=falling_fact_ratio(m,m+l,i)
	fact4=falling_fact_ratio(l,m+l-i,j)
	return fact1*fact2*fact2*fact3*fact4



def do_poc_ratio(a,b,i,j):
	return rising_fact_ratio(1.0*a,1.0*(a+b),i)*rising_fact_ratio(1.0*b,1.0*(a+b+i),j)





def get_sample_sets(nA,nB):
	popC_poly=[]
	popC_maybe_poly=[]
	for mA in range(nA+1):
		for mB in range(nB+1):
			if (mA>0 and mA<nA) and (mB>0 and mB<nB):
				popC_poly.append((mA,nA-mA,mB,nB-mB))
			elif (mA>0 and mA<nA):
				if mB==nB:
					popC_poly.append((mA,nA-mA,nB,0))
				else:
					popC_maybe_poly.append((mA,nA-mA,0,nB))
			elif (mB>0 and mB<nB):
				if mA==nA:
					popC_poly.append((nA,0,mB,nB-mB))
				else:
					popC_maybe_poly.append((0,nA,mB,nB-mB))
			else: #EITHER ONLY DERIVED OR ONLY ANCESTRAL
				popC_maybe_poly.append((mA,nA-mA,mB,nB-mB))
	return [popC_poly,popC_maybe_poly]



def sample_prob_poly(m_A,l_A,m_B,l_B,a,b,tauA,tauB):

	[the_prec,tau_den]=get_prec_den()
	tauA_num=get_tau_num(tauA,tau_den)
	tauB_num=get_tau_num(tauB,tau_den)

	tauA_D=get_tav_dict(m_A+l_A,tauA_num,tau_den,the_prec)
	tauB_D=get_tav_dict(m_B+l_B,tauB_num,tau_den,the_prec)

	h_star_A_dict={}
	for i in range(m_A+1):
		for j in range(l_A+1):
			h_star_A_dict.update({(i,j):h_star(m_A,l_A,i,j)})

	h_star_B_dict={}
	for i in range(m_B+1):
		for j in range(l_B+1):
			h_star_B_dict.update({(i,j):h_star(m_B,l_B,i,j)})

	tot_prob=0.0
	for i_A in range(m_A+1):
		for i_B in range(m_B+1):
			for j_A in range(l_A+1):
				for j_B in range(l_B+1):
					tot_prob+=h_star_A_dict[(i_A,j_A)]*tauA_D['vals'][i_A+j_A]*h_star_B_dict[(i_B,j_B)]*tauB_D['vals'][i_B+j_B]*do_poc_ratio(a,b,i_A+i_B,j_A+j_B)

	return tot_prob





def get_prob_dict_poly(nA,nB,a,b,tauA,tauB):
	[poly,maybe_poly]=get_sample_sets(nA,nB)
	res={}
	totProb=0.0
	for x in poly:
		(mA,lA,mB,lB)=x
		a_val=sample_prob_poly(mA,lA,mB,lB,a,b,tauA,tauB)
		res.update({(mA,lA,mB,lB):a_val})
		totProb+=a_val
	if totProb>0:
		for x in poly:
			res[x]=res[x]/totProb
	return res



def f(x,nA,nB,confs,M):
	[tauA,tauB,a,b]=x
	a_d=get_prob_dict_poly(nA,nB,a,b,tauA,tauB)
	bVal=0.0
	for i in range(len(confs)):
		m=M[i]
		if m>0:
			a_conf=confs[i]
			bVal+=m*np.log(a_d[a_conf])

	#print x,bVal
	return bVal


def f_alt(x):
	[tauA,tauB,a,b]=x
	a_d=get_prob_dict_poly(nA,nB,a,b,tauA,tauB)
	bVal=0.0
	for i in range(len(configurations)):
		m=M[i]
		if m>0:
			a_conf=configurations[i]
			bVal+=m*np.log(a_d[a_conf])

	#print x,bVal
	return -bVal

def find_optimum(nA_test,nB_test,confs_test,M_test):
	global nA, nB, configurations, M
	nA = nA_test
	nB = nB_test
	configurations = confs_test
	M = M_test
	#x_sol = scipy.optimize.basinhopping(f_alt, [0.2, 0.1, 0.005, 1.2])
	bounds =[(0.0000001,100000000000),(0.0000001,10000000000),(0.000000001,1),(0.000000001,1000000)]
	x0 = [0.2, 0.2, 0.005, 1.2]
	x_sol = scipy.optimize.minimize(f_alt, [0.2, 0.8, 0.0005, 1], bounds=bounds)
	#x_sol = scipy.optimize.basinhopping(f_alt, bounds)
	minimizer_kwargs = {"bounds":bounds}
	#x_sol = scipy.optimize.basinhopping(f_alt, x0, minimizer_kwargs=minimizer_kwargs)
	return x_sol
	#return


#global nA, nB, configurations, M
configurations = get_sample_sets(4,4)
nA=4
nB=4
M = [1915, 1123, 714, 351, 1170, 1006, 762, 528, 698, 667, 857, 719, 236, 363, 635]

#res = find_optimum(nA,nB,configurations[0],M)
#res = f([0.1,0.1,0.005,1.2],nA,nB,configs[0],M)
#print(res.x)

# def pop_size(x):
    # return 1/(2*10000)
#
# inter = integrate.quad(pop_size, 0, 2000)
# print(inter)
