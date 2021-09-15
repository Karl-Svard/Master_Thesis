import TT_new_func
import math
import numpy as np

def get_cond_estimates(in_tuple):
    (m10,m01,m20,m02,m11,m21,m12,m00)=in_tuple
    m_tot=1.0*sum(in_tuple)

    if (m10+m20+m21+m11)==0:
        alfa1='NaN'
    else:
        alfa1=1.0*(m10+m12+m11)/(m10+2.0*m20+m21+0.5*m11)

    if (m01+m02+m12+m11)==0:
        alfa2='NaN'
    else:
        alfa2=1.0*(m01+m21+m11)/(m01+2.0*m02+m12+0.5*m11)

    if (m10+m11==0) or (m21+m11==0):
        test1 = 'NaN'
    else:
        test1=(2.0*m10+m11)/(2.0*m01+m11)-(2.0*m12+m11)/(2.0*m21+m11)

    test2=((m10-m01)+2.0*(m20-m02)+(m21-m12))/m_tot
    return [alfa1,alfa2,test1,test2]


def estimate_ancestral_freq_spectra(alfa1,alfa2,p11,p21,p12):
    if alfa1=='NaN' or alfa2=='NaN':
        a_1 = 'NaN'
        a_2 = 'NaN'
        b_1 = 'NaN'
        b_2 = 'NaN'
    else:
        local_s=0.25*p11/(alfa1*alfa2)
        local_r1=0.25*(p11+2.0*p21)/(alfa2)
        local_r2=0.25*(p11+2.0*p12)/(alfa1)

        temp=TT_new_func.find_roots(local_r1, local_s)
        a_1=temp[0]
        b_1=temp[1]

        temp=TT_new_func.find_roots(local_r2, local_s)
        a_2=temp[0]
        b_2=temp[1]

    return [(a_1,b_1),(a_2,b_2)]


def estimate_ancestral_freq_spectra_bound(alfa1,alfa2,p11,p21,p12):
    if alfa1=='NaN' or alfa2=='NaN':
        a_1 = 'NaN'
        a_2 = 'NaN'
        b_1 = 'NaN'
        b_2 = 'NaN'
    else:
        local_s=0.25*p11/(alfa1*alfa2)
        local_r1=0.25*(p11+2.0*p21)/(alfa2)
        local_r2=0.25*(p11+2.0*p12)/(alfa1)

        temp=TT_new_func.find_roots_bounded(local_r1, local_s)
        a_1=temp[0]
        b_1=temp[1]

        temp=TT_new_func.find_roots_bounded(local_r2, local_s)
        a_2=temp[0]
        b_2=temp[1]

    return [(a_1,b_1),(a_2,b_2)]


def estimate_ancestral_freq_spectra_conv(alfa1,alfa2,p11,p21,p12):
    if alfa1=='NaN' or alfa2=='NaN':
        a_1 = 'NaN'
        a_2 = 'NaN'
        b_1 = 'NaN'
        b_2 = 'NaN'
    else:
        local_s=0.25*p11/(alfa1*alfa2)
        local_r1=0.25*(p11+2.0*p21)/(alfa2)
        local_r2=0.25*(p11+2.0*p12)/(alfa1)

        temp=TT_new_func.find_roots(local_r1, local_s)
        temp=TT_new_func.root_converter(temp)
        a_1=temp[0]
        b_1=temp[1]

        temp=TT_new_func.find_roots(local_r2, local_s)
        temp=TT_new_func.root_converter(temp)
        a_2=temp[0]
        b_2=temp[1]

    return [(a_1,b_1),(a_2,b_2)]


def estimate_param(in_tuple,cond_in_tuple,type=''):
    [alfa1,alfa2,test1,test2]=get_cond_estimates(cond_in_tuple)
    (m10,m01,m20,m02,m11,m21,m12,m00)=in_tuple
    m_tot=1.0*sum(in_tuple)
    if type != '':
        estimate_function = globals()['estimate_ancestral_freq_spectra_'+type]
    else:
        estimate_function = globals()['estimate_ancestral_freq_spectra']
    if (alfa1>0.0) and (alfa1<1.0) and (alfa2>0.0) and (alfa2<1.0) and alfa1!='NaN' and alfa2!='NaN':
        [(a_1,b_1),(a_2,b_2)]=estimate_function(alfa1,alfa2,1.0*m11/m_tot,1.0*m21/m_tot,1.0*m12/m_tot)
        if np.isnan(a_1) or np.isnan(a_2):
            a_1='NaN'
            b_1='NaN'
            a_2='NaN'
            b_2='NaN'
            T1_1='NaN'
            T1_2='NaN'
            T2_1='NaN'
            T2_2='NaN'
            V1_1='NaN'
            V1_2='NaN'
            V2_1='NaN'
            V2_2='NaN'
        else:
            x_1=(a_1*b_1)/((a_1+b_1)*(a_1+b_1+1.0))
            x_2=(a_2*b_2)/((a_2+b_2)*(a_2+b_2+1.0))
            y_1=1.0*m10/(2.0*m_tot)+1.0*m11/(4.0*m_tot)+1.0*m20/(1.0*m_tot)+1.0*m21/(2.0*m_tot)
            y_2=1.0*m01/(2.0*m_tot)+1.0*m11/(4.0*m_tot)+1.0*m02/(1.0*m_tot)+1.0*m12/(2.0*m_tot)
            T1_1=y_1-x_1
            T1_2=y_1-x_2
            T2_1=y_2-x_1
            T2_2=y_2-x_2
            V1_1=(2.0*m10-4.0*(alfa1/(1.0-alfa1))*m20-(alfa1*(1.0-alfa1))*((2.0*m21+1.0*m11)+(2.0/(1.0-alfa1))*(m12+m11)))/(4.0*m_tot)
            V2_1=(2.0*m01-4.0*(alfa2/(1.0-alfa2))*m02-(alfa2*(1.0-alfa2))*((2.0*m12+1.0*m11)+(2.0/(1.0-alfa2))*(m21+m11)))/(4.0*m_tot)
            V1_2=(2.0*m10-4.0*(alfa1/(1.0-alfa1))*m20+(1.0/(1.0-alfa1))*m11+(alfa1*(1.0-alfa2)/(alfa2*(1.0-alfa1)))*((2.0*m21+1.0*m11)))/(4.0*m_tot)
            V2_2=(2.0*m01-4.0*(alfa2/(1.0-alfa2))*m02+(1.0/(1.0-alfa2))*m11+(alfa2*(1.0-alfa1)/(alfa1*(1.0-alfa2)))*((2.0*m12+1.0*m11)))/(4.0*m_tot)
    else:
        a_1='NaN'
        b_1='NaN'
        a_2='NaN'
        b_2='NaN'
        T1_1='NaN'
        T1_2='NaN'
        T2_1='NaN'
        T2_2='NaN'
        V1_1='NaN'
        V1_2='NaN'
        V2_1='NaN'
        V2_2='NaN'
    return [alfa1,alfa2,test1,test2,a_1,b_1,a_2,b_2,T1_1,T1_2,T2_1,T2_2,V1_1,V1_2,V2_1,V2_2]


def estimate_param_classic(in_tuple,cond_in_tuple):
	[alfa1,alfa2,test1,test2]=get_cond_estimates(cond_in_tuple)
	(m10,m01,m20,m02,m11,m21,m12,m00)=in_tuple
	m_tot=1.0*sum(in_tuple)
	if alfa1*alfa2>0 and alfa1!='NaN' and alfa2!='NaN':
		y=(9.0*m11)/(m_tot*2.0*alfa1*alfa2)
		tau2_1=(3.0* (2.0*alfa1*m21-(1.0-alfa1)*m11) ) /(m_tot*2.0*alfa1*alfa2)
		tau2_2=(3.0* (2.0*alfa2*m12-(1.0-alfa2)*m11) ) /(m_tot*2.0*alfa1*alfa2)
		tau3_1=((5.0-2.0*alfa1)*m11-4.0*alfa1*m21 ) /(m_tot*2.0*alfa1*alfa2)
		tau3_2=((5.0-2.0*alfa2)*m11-4.0*alfa2*m12 ) /(m_tot*2.0*alfa1*alfa2)
		B1=(0.5*m10+m20+0.5*m21+0.25*m11-(5.0*m11/(4.0*alfa1*alfa2)) )/m_tot
		B2=(0.5*m01+m02+0.5*m12+0.25*m11-(5.0*m11/(4.0*alfa1*alfa2)) )/m_tot
		#T1=B1-(y/18.0)
		#T2=B2-(y/18.0)
		#T1=B1-0.25*(tau3_1+tau3_2)
		#T2=B2-0.25*(tau3_1+tau3_2)
		#mean_B=0.5*(B1+B2)
		mean_tau_3=0.5*(tau3_1+tau3_2)
		mean_tau_2=0.5*(tau2_1+tau2_2)
		T1=B1-0.5*mean_tau_3
		T2=B2-0.5*mean_tau_3
		J1=B1-(3.0*mean_tau_3)*(3.0*mean_tau_3)/(6.0*mean_tau_2)
		J2=B2-(3.0*mean_tau_3)*(3.0*mean_tau_3)/(6.0*mean_tau_2)
	else:
		y='NaN'
		tau2_1='NaN'
		tau2_2='NaN'
		tau3_1='NaN'
		tau3_2='NaN'
		B1='NaN'
		B2='NaN'
		T1='NaN'
		T2='NaN'
		J1='NaN'
		J2='NaN'

	if alfa1!='NaN' and alfa2!='NaN':
		U1=(0.5*(1.0-alfa1)*m10-alfa1*m20+0.5*(1.0-alfa2)*m12 +0.25*(2.0-alfa2)*m11 )/m_tot
		U2=(0.5*(1.0-alfa2)*m01-alfa2*m02+0.5*(1.0-alfa1)*m21 +0.25*(2.0-alfa1)*m11 )/m_tot
	else:
		U1='Nan'
		U2='NaN'

	if alfa1<1 and alfa1!='NaN':
		#V1=(0.5*m10-(alfa1/(1.0-alfa1))*m20+0.5*((1.0-alfa2)/(1.0-alfa1))*m12 +0.25*((2.0-alfa2)/(1.0-alfa1))*m11 )/m_tot
		V1=(0.5*m10-  (alfa1/(1.0-alfa1))*m20  +0.5*m12/(1.0-alfa1)  )/m_tot
	else:
		V1='NaN'

	if alfa2<1 and alfa2!='NaN':
		#V2=(0.5*m01-(alfa2/(1.0-alfa2))*m02+0.5*((1.0-alfa1)/(1.0-alfa2))*m21 +0.25*((2.0-alfa1)/(1.0-alfa2))*m11 )/m_tot
		V2=(0.5*m01-(alfa2/(1.0-alfa2))*m02+0.5*m21/(1.0-alfa2))/m_tot
	else:
		V2='NaN'

	if tau2_1 != 'NaN' and tau2_2 != 'NaN':
		tau_test=y-1.5*(tau2_1+tau2_2)
	else:
		tau_test = 'NaN'

	return [alfa1,alfa2,test1,test2,y,tau2_1,tau2_2,tau3_1,tau3_2,B1,B2,U1,U2,V1,V2,tau_test,T1,T2,J1,J2]


def estimate_param_TT(in_tuple):
    (n1,n2,n3,n4,n5,n6,n7,n0)=in_tuple
    n_tot=1.0*sum(in_tuple)

    if (n5+n6)==0:
        alfa1='NaN'
    else:
        alfa1=2.0*n5/(n5+2.0*n6)
    if (n5+n7)==0:
        alfa2='NaN'
    else:
        alfa2=2.0*n5/(n5+2.0*n7)

    mu_t1_t2_diff=1.0*(n3-n4+0.5*(n1+n6-n2-n7))/n_tot

    if n5<1:
        thetaA='NaN'
        mu_t1='NaN'
        mu_t2='NaN'
    else:
        thetaA=(3.0/n_tot)*(n5+2.0*n6)*(n5+2.0*n7)/(8.0*n5)
        mu_t1_p1=1.0*n3+0.5*n1
        mu_t1_p2=(n5+2.0*n6)*(n5+6.0*n7)/(8.0*n5)
        mu_t1=(mu_t1_p1-mu_t1_p2)/n_tot

        mu_t2_p1=1.0*n4+0.5*n2
        mu_t2_p2=(n5+6.0*n6)*(n5+2.0*n7)/(8.0*n5)
        mu_t2=(mu_t2_p1-mu_t2_p2)/n_tot

    if n5<2*n6:
	#drift1=-1.0*math.log(alfa1)
	#theta1=mu_t1/drift1
        the_nom=2.0*n6*(n1+n7)-1.0*n5*(n1+4.0*n3-n7)
        the_den=2.0*(2.0*n6-1.0*n5)*n_tot
        mu_nu1=the_nom/the_den
        the_nom=(2.0*n6+n5)*(n5*(8.0*n3+2.0*n7+n5)-2.0*n6*(6.0*n7+n5))
        the_den=(2.0*n6-n5)*(n5*(4.0*n1+8.0*n3-6.0*n7-n5)-2.0*n6*(6.0*n7+n5))
        
        if the_den>0:
        	W1ratio=the_nom/the_den
        else:
        	W1ratio='NaN'

        if n5>0:
    	        logpt=1.0/(math.log(2.0*n5)-math.log(2.0*n6+n5))
    	        pt1=(2.0*n7+n1)*(2.0*n6+n5)/(2.0*n6-n5)
    	        pt2=(2.0*n3-n1)/logpt
    	        pt3=(2.0*n6+n5)*(6.0*n7+n5)/(4.0*n5)
    	        pt4=4.0*n5/(2.0*n6-n5)
    	        D1=(pt1+pt2-(pt3*(pt4+logpt)))/(2.0*n_tot)
        else:
                D1='NaN'

        if alfa1 != 'NaN' and mu_t1 != 'NaN':
                drift1=-1.0*math.log(alfa1)
                theta1=mu_t1/drift1
        else:
                drift1='NaN'
                theta1='NaN'

    else:
        drift1='NaN'
        theta1='NaN'
        mu_nu1='NaN'
        W1ratio='NaN'
        D1='NaN'

    if n5<2*n7:
	#drift2=-1.0*math.log(alfa2)
	#theta2=mu_t2/drift2
        the_nom=2.0*n7*(n2+n6)-1.0*n5*(n2+4.0*n4-n6)
        the_den=2.0*(2.0*n7-1.0*n5)*n_tot
        mu_nu2=the_nom/the_den

        the_nom=(2.0*n7+n5)*(n5*(8.0*n4+2.0*n6+n5)-2.0*n7*(6.0*n6+n5))
        the_den=(2.0*n7-n5)*(n5*(4.0*n2+8.0*n4-6.0*n6-n5)-2.0*n7*(6.0*n6+n5))
        if the_den>0:
        	W2ratio=the_nom/the_den
        else:
        	W2ratio='NaN'

        if n5>0:
        	logpt=1.0/(math.log(2.0*n5)-math.log(2.0*n7+n5))
        	pt1=(2.0*n6+n2)*(2.0*n7+n5)/(2.0*n7-n5)
        	pt2=(2.0*n4-n2)/logpt
        	pt3=(6.0*n6+n5)*(2.0*n7+n5)/(4.0*n5)
        	pt4=4.0*n5/(2.0*n7-n5)
        	D2=(pt1+pt2-(pt3*(pt4+logpt)))/(2.0*n_tot)
        else:
            D2='NaN'

        if alfa2 != 'NaN' and mu_t2 != 'NaN':
        	drift2=-1.0*math.log(alfa2)
        	theta2=mu_t2/drift2
        else:
                drift2='NaN'
                theta2='NaN'

    else:
    	drift2='NaN'
    	theta2='NaN'
    	mu_nu2='NaN'
    	W2ratio='NaN'
    	D2='NaN'
##################METHOD Schlebusch et al 2012 ############
    if (n3+n5+n6+n7)>0 and thetaA != 'NaN':
        Ccount1=1.0*n3+0.5*n6
        D1orD2count1=0.5*n5+1.0*n7
        P1=-math.log((3.0/2.0)*D1orD2count1/(D1orD2count1+Ccount1))
        P1_time=P1*thetaA
    else:
        P1 = 'NaN'
        P1_time = 'NaN'

    if (n4+n5+n6+n7)>0 and thetaA != 'NaN':
    	Ccount2=1.0*n4+0.5*n7
    	D1orD2count2=0.5*n5+1.0*n6
    	P2=-math.log((3.0/2.0)*D1orD2count2/(D1orD2count2+Ccount2))
    	P2_time=P2*thetaA
    else:
        P2 = 'NaN'
        P2_time = 'NaN'

###################METHOD Schlebusch et al 2012 ############
###################Fst############
    if n_tot!=n0:
        Fst=(2.0*n3+2.0*n4-1.0*n5)/(1.0*n1+1.0*n2+2.0*n3+2.0*n4+1.0*n5+1.0*n6+1.0*n7)
    else:
        Fst='NaN'
###################Fst############
    return [alfa1,alfa2,thetaA,mu_t1,mu_t2,mu_nu1,mu_nu2,mu_t1_t2_diff,drift1,drift2,theta1,theta2,W1ratio,W2ratio,D1,D2,P1,P2,P1_time,P2_time,Fst]
