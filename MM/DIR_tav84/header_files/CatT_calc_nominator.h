
#ifndef CATT_NOMINATOR
#define CATT_NOMINATOR

int function_nominator(mpfr_t mpfr_BACK, ListNode* SENT_tav_vect, bool prior_uniform,bool U_is_derived, unsigned long int m_2,unsigned long int l_2,unsigned long int m_3,unsigned long int l_3,unsigned long int m_A,unsigned long int l_A,unsigned long int c_nom,unsigned long int c_den,unsigned long int t_ac_nom,unsigned long int t_ac_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC);

#endif


int function_nominator(mpfr_t mpfr_BACK, ListNode* SENT_tav_vect, bool prior_uniform,bool U_is_derived, unsigned long int m_2,unsigned long int l_2,unsigned long int m_3,unsigned long int l_3,unsigned long int m_A,unsigned long int l_A,unsigned long int c_nom,unsigned long int c_den,unsigned long int t_ac_nom,unsigned long int t_ac_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC)
{

  mpfr_t mpfr_A,mpfr_B,mpfr_B2;
  mpq_t mpq_C,mpq_D,mpq_E,mpq_F;
  mpz_t mpz_G;

  unsigned long int n;
  unsigned long int alpha_2,beta_2,alpha_3,beta_3,i_2,j_2,i_3,j_3;
  unsigned long int tav_alpha,tav_j;
  unsigned long int temp_ui_1,temp_ui_2;

  ListNode* Mpfr_tav_vect=SENT_tav_vect;
  ListNode* at_node=NULL;
  ListNode* a_new_node=NULL;
  ListNode* temp_node=NULL;
  mpfr_t* an_mpfr_ptr=NULL;



//INITIALISE MPFR AND GMP VARIABLES

  mpfr_init2(mpfr_A,MY_PREC);
  mpfr_init2(mpfr_B,MY_PREC);
  mpfr_init2(mpfr_B2,MY_PREC);
  mpq_init(mpq_C);
  mpq_init(mpq_D);
  mpq_init(mpq_E);
  mpq_init(mpq_F);
  mpz_init(mpz_G);

  n=m_2+l_2;

 //HERE ARE CALCULATIONS FOR THE DENOMINATOR WITH STUCTURE:
  /*
    A=0  (MPFR)
    for (alpha_2=0;alpha_2<=m_2;alpha_2++)
      {
	for (beta_2=0;beta_2<=l_2;beta_2++)
	  {
            for (alpha_3=0;alpha_3<=m_3;alpha_3++)
              {
	        for (beta_3=0;beta_3<=l_3;beta_3++)
	          {
	            C=0 (GMP)
	            for (i_2=0;i_2<=alpha_2;i_2++)
	              {
                        for (j_2=0;j_2<=beta_2;j_2++)
                          {
         	            for (i_3=0;i_3<=alpha_3;i_3++)
	                      {
                                for (j_3=0;j_3<=beta_3;j_3++)
                                  {
		                    C=C+func_my_Bin(alpha_2,i_2,c)*func_my_Bin(beta_2,j_2,c)*func_my_Bin(alpha_3,i_3,c)*func_my_Bin(beta_3,j_3,c)
						*func_special_Beta(prior_uniform,U_is_derived,m_A+i_2,l_a+j_2)*func_special_Beta(prior_uniform,U_is_derived,m_A+i_3,l_a+j_3)
							*func_Beta(alpha_2+alpha_3-i_2-i_3+1,beta_2+beta_3-j_2-j_3+1)
                                  }//j_3-loop
                              }//i_3-loop
                          }//j_2-loop
                      }//i_2-loop
	            B=C*func_n_choose_k(alpha_2+beta_2,alpha_2)*func_n_choose_k(alpha_3+beta_3,alpha_3)
			*func_H(alpha_2,beta_2,m_2,l_2)*func_H(alpha_3,beta_3,m_3,l_3)  (GMP)
	            A=A+B*func_G(n,alpha_2+beta_2,t)*func_G(n,alpha_3+beta_3) (MPFR)
                  }//beta_3-loop
              }//alpha_3-loop
          }//beta_2-loop
      }//alpha_2-loop
    C=func_special_Beta(prior_uniform,U_is_derived,m_A,l_A)*func_special_Beta(prior_uniform,U_is_derived,m_A,l_A)
    A=A/C
  */


  mpfr_set_ui(mpfr_A,0,ROUND_MODE);
  for (alpha_2=0;alpha_2<=m_2;alpha_2++)
    {
      for (beta_2=0;beta_2<=l_2;beta_2++)
        {
          for (alpha_3=0;alpha_3<=m_3;alpha_3++)
            {
              for (beta_3=0;beta_3<=l_3;beta_3++)
                {
                  mpq_set_ui(mpq_D,0,1);
                  for (i_2=0;i_2<=alpha_2;i_2++)
	            {
                      for (j_2=0;j_2<=beta_2;j_2++)
                        {
                  	  for (i_3=0;i_3<=alpha_3;i_3++)
	            	    {
                      	      for (j_3=0;j_3<=beta_3;j_3++)
                        	{
		  		  function_myBin(mpq_E,alpha_2,i_2,c_nom,c_den);
		  		  function_myBin(mpq_F,beta_2,j_2,c_nom,c_den);
		  	   	  mpq_mul(mpq_E,mpq_E,mpq_F);
		  		  mpq_canonicalize(mpq_E);
		  		  function_myBin(mpq_F,alpha_3,i_3,c_nom,c_den);
		  		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  		  mpq_canonicalize(mpq_E);
		  		  function_myBin(mpq_F,beta_3,j_3,c_nom,c_den);
		  		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  		  mpq_canonicalize(mpq_E);
		  		  function_special_Beta(mpq_F,prior_uniform,U_is_derived,m_A+i_2,l_A+j_2);
		  		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  		  mpq_canonicalize(mpq_E);
		  		  function_special_Beta(mpq_F,prior_uniform,U_is_derived,m_A+i_3,l_A+j_3);
		  		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  		  mpq_canonicalize(mpq_E);
		  		  temp_ui_1=alpha_2+alpha_3-i_2-i_3+1;
		  		  temp_ui_2=beta_2+beta_3-j_2-j_3+1;
		  		  function_Beta(mpq_F,temp_ui_1,temp_ui_2);
		  		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  		  mpq_canonicalize(mpq_E);
	          		  mpq_add(mpq_D,mpq_D,mpq_E);
		  		  mpq_canonicalize(mpq_D);
                		}//j_3-loop
	    		    }//i_3-loop
                	}//j_2-loop
	            }//i_2-loop
	          temp_ui_1=alpha_2+beta_2;
	          mpz_bin_uiui(mpz_G,temp_ui_1,alpha_2);
	          mpq_set_z(mpq_C,mpz_G);
	  	  mpq_canonicalize(mpq_C);
	  	  function_H(mpq_E,alpha_2,beta_2,m_2,l_2);
	  	  mpq_mul(mpq_C,mpq_C,mpq_E);
	  	  mpq_canonicalize(mpq_C);

	          temp_ui_1=alpha_3+beta_3;
	          mpz_bin_uiui(mpz_G,temp_ui_1,alpha_3);
	          mpq_set_z(mpq_F,mpz_G);
	  	  function_H(mpq_E,alpha_3,beta_3,m_3,l_3);
		  mpq_mul(mpq_E,mpq_E,mpq_F);
	  	  mpq_canonicalize(mpq_E);
		  mpq_mul(mpq_C,mpq_C,mpq_E);
	  	  mpq_canonicalize(mpq_C);
	  	  mpq_mul(mpq_C,mpq_C,mpq_D);
	  	  mpq_canonicalize(mpq_C);

          	  temp_node=Mpfr_tav_vect->GetBetaAlphaNode(Mpfr_tav_vect,n,alpha_2+beta_2);
          	  an_mpfr_ptr=temp_node->GetMpfr();
	  	  mpfr_set(mpfr_B,*an_mpfr_ptr,ROUND_MODE);
          	  temp_node=Mpfr_tav_vect->GetBetaAlphaNode(Mpfr_tav_vect,n,alpha_3+beta_3);
          	  an_mpfr_ptr=temp_node->GetMpfr();
	  	  mpfr_set(mpfr_B2,*an_mpfr_ptr,ROUND_MODE);
	  	  mpfr_mul(mpfr_B,mpfr_B,mpfr_B2,ROUND_MODE);
	  	  mpfr_mul_q(mpfr_B,mpfr_B,mpq_C,ROUND_MODE);
	  	  temp_node=NULL;
	  	  an_mpfr_ptr=NULL;

	  	  mpfr_add(mpfr_A,mpfr_A,mpfr_B,ROUND_MODE);
		} //beta_3-loop
    	    }//alpha_3-loop
	} //beta_2-loop
    }//alpha_2-loop

  function_special_Beta(mpq_E,prior_uniform,U_is_derived,m_A,l_A);
  function_special_Beta(mpq_F,prior_uniform,U_is_derived,m_A,l_A);
  mpq_mul(mpq_E,mpq_E,mpq_F);
  mpq_canonicalize(mpq_E);
  mpfr_div_q(mpfr_A,mpfr_A,mpq_E,ROUND_MODE);
  mpfr_swap(mpfr_BACK,mpfr_A);
  mpq_clear(mpq_F);
  mpq_clear(mpq_E);
  mpq_clear(mpq_D);
  mpq_clear(mpq_C);
  mpfr_clear(mpfr_B2);
  mpfr_clear(mpfr_B);
  mpfr_clear(mpfr_A);
  mpfr_free_cache();
  return 0;
}

