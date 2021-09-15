#ifndef MAIN_PROBABILITIES
#define MAIN_PROBABILITIES


int function_Beta(mpq_t gmp_dummy,unsigned long int n,unsigned long int k);
int function_special_Beta(mpq_t gmp_dummy,bool prior_uniform, bool U_is_derived, unsigned long int n,unsigned long int k);
int function_myBin(mpq_t gmp_dummy,unsigned long int n,unsigned long int k,unsigned long int a_nom,unsigned long int a_den);
//int function_C(mpq_t gmp_dummy,unsigned long int n,unsigned long int k,unsigned long int c_nom,unsigned long int c_den);
//function_C and function_Bin are identical and return (n over k)c^(k)(1-c)^(n-k)
int check_which_case_H_func(unsigned long int j,unsigned long int i,unsigned long int m,unsigned long int k);
int function_H(mpq_t gmp_dummy,unsigned long int j,unsigned long int i,unsigned long int m,unsigned long int k);
int check_which_case_G_func(unsigned long int m,unsigned long int i);
int function_G(mpfr_t mpfr_dummy,unsigned long int m,unsigned long int i,unsigned long int t_nom,unsigned long int t_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC);

#endif

int function_Beta(mpq_t gmp_dummy,unsigned long int n,unsigned long int k)
{
  mpz_t gmp_Local_A, gmp_Local_B;
  unsigned long int n_minus_1,k_minus_1,n_plus_k_minus_1;

  n_minus_1=n-1;
  k_minus_1=k-1;
  n_plus_k_minus_1=n+k-1;

  mpz_init(gmp_Local_A);
  mpz_init(gmp_Local_B);

  mpz_fac_ui(gmp_Local_A,n_minus_1);
  mpz_fac_ui(gmp_Local_B,k_minus_1);
  mpz_mul(gmp_Local_A,gmp_Local_A,gmp_Local_B);
  mpz_fac_ui(gmp_Local_B,n_plus_k_minus_1);

  mpq_set_num(gmp_dummy,gmp_Local_A);
  mpq_set_den(gmp_dummy,gmp_Local_B);
  mpq_canonicalize(gmp_dummy);
  mpz_clear(gmp_Local_A);
  mpz_clear(gmp_Local_B);
  mpfr_free_cache();
  return 0;
}

int function_special_Beta(mpq_t gmp_dummy,bool prior_uniform, bool U_is_derived, unsigned long int n,unsigned long int k)
{
  unsigned long int n_plus_1,k_plus_1;
  n_plus_1=n+1;
  k_plus_1=k+1;

  if (prior_uniform==true)
    {
      function_Beta(gmp_dummy,n_plus_1,k_plus_1);
    }
  else
    {
      if (U_is_derived==true)
        {
          function_Beta(gmp_dummy,n,k_plus_1);
        }
      else
        {
          function_Beta(gmp_dummy,n_plus_1,k);
        }
    }
  return 0;
}


int function_myBin(mpq_t gmp_dummy,unsigned long int n,unsigned long int k,unsigned long int a_nom,unsigned long int a_den)
{
  mpz_t gmp_Local_A, gmp_Local_B;
  unsigned long int a_diff,n_minus_k;

  a_diff=a_den-a_nom;
  n_minus_k=n-k;
  mpz_init(gmp_Local_A);
  mpz_init(gmp_Local_B);

  mpz_ui_pow_ui(gmp_Local_A,a_nom,k);
  mpz_ui_pow_ui(gmp_Local_B,a_diff,n_minus_k);
  mpz_mul(gmp_Local_A,gmp_Local_A,gmp_Local_B);
  mpz_bin_uiui(gmp_Local_B,n,k);
  mpz_mul(gmp_Local_A,gmp_Local_A,gmp_Local_B);
  mpz_ui_pow_ui(gmp_Local_B,a_den,n);

  mpq_set_num(gmp_dummy,gmp_Local_A);
  mpq_set_den(gmp_dummy,gmp_Local_B);
  mpq_canonicalize(gmp_dummy);
  mpz_clear(gmp_Local_A);
  mpz_clear(gmp_Local_B);
  mpfr_free_cache();
  return 0;
}


/* THIS DOESN'T WORK FOR SOME REASON
int function_C(mpq_t gmp_dummy,unsigned long int n,unsigned long int k,unsigned long int c_nom,unsigned long int c_den)
{
  function_Bin(mpq_t gmp_dummy,unsigned long int n,unsigned long int k,unsigned long int c_nom,unsigned long int c_den);
  return 0;
}
*/

int check_which_case_H_func(unsigned long int j,unsigned long int i,unsigned long int m,unsigned long int k)
{
  int ret_val;
  if (((j>0)&&(j<=m))&&((i>0)&&(i<=k)))
    {
      ret_val=2;
    }
  else
    {
      if (((j==0)&&(m==0))&&((i>0)&&(i<=k)))
	{
	  ret_val=1;
	}
      else
	{
	  if (((j>0)&&(j<=m))&&((i==0)&&(k==0)))
	    {
	      ret_val=1;
	    }
	  else
	    {
	      if (((j==0)&&(m==0))&&((i==0)&&(k==0)))
		{
		  ret_val=1;
		}
	      else
		{
		  ret_val=0;
		}
	    }
	}
    }
  return ret_val;
}

int function_H(mpq_t gmp_dummy,unsigned long int j,unsigned long int i,unsigned long int m,unsigned long int k)
{
  int case_indicator,ret_val=1;
  mpz_t gmp_Local_A, gmp_Local_B;
  mpz_init(gmp_Local_A);
  mpz_init(gmp_Local_B);
  case_indicator=check_which_case_H_func(j,i,m,k);
  if (case_indicator==0)
    {
      mpz_set_ui(gmp_Local_A,0);
      mpz_set_ui(gmp_Local_B,1);
      mpq_set_num(gmp_dummy,gmp_Local_A);
      mpq_set_den(gmp_dummy,gmp_Local_B);
      mpq_canonicalize(gmp_dummy);
    }
  else
    {
       if (case_indicator==1)
	 {
	   mpz_set_ui(gmp_Local_A,1);
	   mpz_set_ui(gmp_Local_B,1);
	   mpq_set_num(gmp_dummy,gmp_Local_A);
	   mpq_set_den(gmp_dummy,gmp_Local_B);
	   mpq_canonicalize(gmp_dummy);
	 }
       else
	 {
	   mpz_bin_uiui(gmp_Local_A,m-1,j-1);
	   mpz_bin_uiui(gmp_Local_B,k-1,i-1);
	   mpz_mul(gmp_Local_A,gmp_Local_A,gmp_Local_B);
	   mpz_bin_uiui(gmp_Local_B,m+k-1,j+i-1);
	   mpq_set_num(gmp_dummy,gmp_Local_A);
	   mpq_set_den(gmp_dummy,gmp_Local_B);
	   mpq_canonicalize(gmp_dummy);
	 }
    }
  mpz_clear(gmp_Local_A);
  mpz_clear(gmp_Local_B);
  mpfr_free_cache();
  return ret_val;
}


int check_which_case_G_func(unsigned long int m,unsigned long int i,unsigned long int t_nom)
{
  int ret_val;
  if (!(t_nom==0))
    {
      if ((0<i)&&(i<=m))
	{
	  ret_val=2;
	}
      else
	{
	  if ((i==0)&&(m==0))
	    {
	      ret_val=1;
	    }
	  else
	    {
	      ret_val=0;
	    }
	}
    }
  else
    {
      if (i==m)
	{
	  ret_val=1;
	}
      else
	{
	  ret_val=0;
	}
    }
  return ret_val;
}


int function_G(mpfr_t mpfr_dummy,unsigned long int m,unsigned long int i,unsigned long int t_nom,unsigned long int t_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC)
{
  int  case_indicator,non_negative,not_larger_than_one,ret_val=1;
  unsigned long int local_r;
  mpz_t gmp_Local_A, gmp_Local_B;
  mpfr_t mpfr_Local;

  case_indicator=check_which_case_G_func(m,i,t_nom);
  if (case_indicator==0)
    {
      mpfr_set_ui(mpfr_dummy,0,ROUND_MODE);
    }
  else
    {
      if (case_indicator==1)
	{
	  mpfr_set_ui(mpfr_dummy,1,ROUND_MODE);
	}
      else
	{
	  mpz_init(gmp_Local_A);
	  mpz_init(gmp_Local_B);
	  mpfr_init2(mpfr_Local,MY_PREC);
	  mpfr_set_ui(mpfr_dummy,0,ROUND_MODE);
	  if (i==1)
	    {
	      for (local_r=2; local_r<=m; local_r++)
		{
		  mpfr_set_ui(mpfr_Local,(local_r*local_r-local_r)*t_nom,ROUND_MODE);
		  mpfr_div_ui(mpfr_Local,mpfr_Local,2*t_den,ROUND_MODE);
		  mpfr_neg(mpfr_Local,mpfr_Local,ROUND_MODE);
		  mpfr_exp(mpfr_Local,mpfr_Local,ROUND_MODE);
		  mpfr_mul_ui(mpfr_Local,mpfr_Local,2*local_r-1,ROUND_MODE);
		  if (local_r%2==1)
		    {
		      mpfr_neg(mpfr_Local,mpfr_Local,ROUND_MODE);
		    }
		  mpz_bin_uiui(gmp_Local_A,m,local_r);
		  mpz_bin_uiui(gmp_Local_B,m+local_r-1,local_r);
		  mpfr_mul_z(mpfr_Local,mpfr_Local,gmp_Local_A,ROUND_MODE);
		  mpfr_div_z(mpfr_Local,mpfr_Local,gmp_Local_B,ROUND_MODE);
		  mpfr_add(mpfr_dummy,mpfr_dummy,mpfr_Local,ROUND_MODE);
		}
	      mpfr_neg(mpfr_dummy,mpfr_dummy,ROUND_MODE);
	      mpfr_add_ui(mpfr_dummy,mpfr_dummy,1,ROUND_MODE);
	    }
	  else
	    {
	      for (local_r=i; local_r<=m; local_r++)
		{
		  mpfr_set_ui(mpfr_Local,(local_r*local_r-local_r)*t_nom,ROUND_MODE);
		  mpfr_div_ui(mpfr_Local,mpfr_Local,2*t_den,ROUND_MODE);
		  mpfr_neg(mpfr_Local,mpfr_Local,ROUND_MODE);
		  mpfr_exp(mpfr_Local,mpfr_Local,ROUND_MODE);
		  mpfr_mul_ui(mpfr_Local,mpfr_Local,2*local_r-1,ROUND_MODE);
		  if ((local_r-i)%2==1)
		    {
		      mpfr_neg(mpfr_Local,mpfr_Local,ROUND_MODE);
		    }
		  mpz_bin_uiui(gmp_Local_A,m,local_r);

		  mpz_bin_uiui(gmp_Local_B,local_r,i);
		  mpz_mul(gmp_Local_A,gmp_Local_A,gmp_Local_B);

		  mpz_bin_uiui(gmp_Local_B,i+local_r-1,local_r);
		  mpz_mul(gmp_Local_A,gmp_Local_A,gmp_Local_B);

		  mpz_bin_uiui(gmp_Local_B,m+local_r-1,local_r);
		  mpz_mul_ui(gmp_Local_B,gmp_Local_B,i+local_r-1);

		  mpfr_mul_z(mpfr_Local,mpfr_Local,gmp_Local_A,ROUND_MODE);
		  mpfr_div_z(mpfr_Local,mpfr_Local,gmp_Local_B,ROUND_MODE);
		  mpfr_add(mpfr_dummy,mpfr_dummy,mpfr_Local,ROUND_MODE);
		}
	    }
	  non_negative=mpfr_cmp_ui(mpfr_dummy,0);
	  not_larger_than_one=-1*mpfr_cmp_ui(mpfr_dummy,1);
	  if ((non_negative<0)||(not_larger_than_one<0))
	    {
	      ret_val=0;
	    }
	  mpfr_clear(mpfr_Local);
	  mpz_clear(gmp_Local_A);
	  mpz_clear(gmp_Local_B);
	  mpfr_free_cache();
	}
    }
  return ret_val;
}

