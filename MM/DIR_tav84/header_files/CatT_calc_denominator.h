
#ifndef CATT_DENOMINATOR
#define CATT_DENOMINATOR

int function_denominator(mpfr_t mpfr_BACK, ListNode* SENT_tav_vect, bool prior_uniform, bool U_is_derived, unsigned long int m_C,unsigned long int l_C, unsigned long int m_A,unsigned long int l_A,unsigned long int c_nom,unsigned long int c_den,unsigned long int t_ac_nom,unsigned long int t_ac_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC);

#endif

int function_denominator(mpfr_t mpfr_BACK, ListNode* SENT_tav_vect, bool prior_uniform, bool U_is_derived, unsigned long int m_C,unsigned long int l_C, unsigned long int m_A,unsigned long int l_A,unsigned long int c_nom,unsigned long int c_den,unsigned long int t_ac_nom,unsigned long int t_ac_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC)
{


  ListNode* Mpfr_tav_vect=SENT_tav_vect;
  ListNode* at_node=NULL;
  ListNode* a_new_node=NULL;
  ListNode* temp_node=NULL;
  mpfr_t* an_mpfr_ptr=NULL;


  mpfr_t mpfr_A,mpfr_B;
  mpq_t mpq_C,mpq_D,mpq_E,mpq_F;
  mpz_t mpz_G;

  unsigned long int n;
  unsigned long int temp_ui_1,temp_ui_2;
  unsigned long int tav_alpha,tav_j,i,j,alpha,beta;



//INITIALISE MPFR AND GMP VARIABLES
  mpfr_init2(mpfr_A,MY_PREC);
  mpfr_init2(mpfr_B,MY_PREC);
  mpq_init(mpq_C);
  mpq_init(mpq_D);
  mpq_init(mpq_E);
  mpq_init(mpq_F);
  mpz_init(mpz_G);


  n=m_C+l_C;

  ///////////////////////////////
 //HERE ARE CALCULATIONS FOR THE DENOMINATOR WITH STUCTURE:
  /*
    A=0  (MPFR)
    for (alpha=0;alpha<=m_C;alpha++)
      {
	for (beta=0;beta<=l_C;beta++)
	  {
            C=0 (GMP)
	    for (i=0;i<=alpha;i++)
	      {
                for (j=0;j<=beta;j++)
                  {
		     C=C+func_my_Bin(alpha,i,c)*func_my_Bin(beta,j,c)*func_special_Beta(prior_uniform,U_is_derived,m_A+i,l_a+j)*func_Beta(alpha-i+1,beta-j+1)
                  }//j-loop
              }//i-loop
	    B=C*func_n_choose_k(alpha+beta,alpha)*func_H(alpha,beta,m_C,l_C) (GMP)
	    A=A+B*func_G(n,alpha+beta,t) (MPFR)
          }//beta-loop
       }//alpha-loop
    A=A/func_special_Beta(m_A,l_A)
  */
  mpfr_set_ui(mpfr_A,0,ROUND_MODE);
  for (alpha=0;alpha<=m_C;alpha++)
    {
      for (beta=0;beta<=l_C;beta++)
        {
	  mpq_set_ui(mpq_D,0,1);
          for (i=0;i<=alpha;i++)
	    {
              for (j=0;j<=beta;j++)
                {
		  temp_ui_1=m_A+i;
		  temp_ui_2=l_A+j;
		  function_special_Beta(mpq_E,prior_uniform,U_is_derived,temp_ui_1,temp_ui_2);
		  function_myBin(mpq_F,alpha,i,c_nom,c_den);
		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  mpq_canonicalize(mpq_E);
		  function_myBin(mpq_F,beta,j,c_nom,c_den);
		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  mpq_canonicalize(mpq_E);
		  temp_ui_1=alpha-i+1;
		  temp_ui_2=beta-j+1;
		  function_Beta(mpq_F,temp_ui_1,temp_ui_2);
		  mpq_mul(mpq_E,mpq_E,mpq_F);
		  mpq_canonicalize(mpq_E);
	          mpq_add(mpq_D,mpq_D,mpq_E);
		  mpq_canonicalize(mpq_D);
                }//j-loop
	    }//i-loop
	  temp_ui_1=alpha+beta;
	  mpz_bin_uiui(mpz_G,temp_ui_1,alpha);
	  mpq_set_z(mpq_C,mpz_G);
	  function_H(mpq_E,alpha,beta,m_C,l_C);
	  mpq_mul(mpq_C,mpq_C,mpq_E);
	  mpq_canonicalize(mpq_C);
	  mpq_mul(mpq_C,mpq_C,mpq_D);
	  mpq_canonicalize(mpq_C);
          temp_node=Mpfr_tav_vect->GetBetaAlphaNode(Mpfr_tav_vect,n,temp_ui_1);
          an_mpfr_ptr=temp_node->GetMpfr();
	  mpfr_set(mpfr_B,*an_mpfr_ptr,ROUND_MODE);
	  temp_node=NULL;
	  an_mpfr_ptr=NULL;
	  mpfr_mul_q(mpfr_B,mpfr_B,mpq_C,ROUND_MODE);
	  mpfr_add(mpfr_A,mpfr_A,mpfr_B,ROUND_MODE);
	} //beta-loop
    }//alpha-loop
  function_special_Beta(mpq_E,prior_uniform,U_is_derived,m_A,l_A);
  mpfr_div_q(mpfr_A,mpfr_A,mpq_E,ROUND_MODE);
  mpfr_swap(mpfr_BACK,mpfr_A);
  mpq_clear(mpq_F);
  mpq_clear(mpq_E);
  mpq_clear(mpq_D);
  mpq_clear(mpq_C);
  mpfr_clear(mpfr_B);
  mpfr_clear(mpfr_A);
  mpfr_free_cache();
  return 0;
}

