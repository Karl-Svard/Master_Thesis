#ifndef TAV_LISTS
#define TAV_LISTS

ListNode* function_Tav_1D(unsigned long int n,unsigned long int t_ac_nom,unsigned long int t_ac_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC);
#endif


ListNode* function_Tav_1D(unsigned long int n,unsigned long int t_ac_nom,unsigned long int t_ac_den, mp_rnd_t ROUND_MODE,mp_prec_t MY_PREC)
{


  ListNode* Mpfr_tav_vect;
  ListNode* at_node=NULL;
  ListNode* a_new_node=NULL;
  ListNode* temp_node=NULL;
  mpfr_t* an_mpfr_ptr=NULL;
  mpfr_t mpfr_A;
  unsigned long int tav_alpha,tav_j;

  mpfr_init2(mpfr_A,MY_PREC);
  Mpfr_tav_vect=new ListNode(n,0,MY_PREC);
  at_node=Mpfr_tav_vect;
  if (n==0)
    {
      mpfr_set_ui(mpfr_A,1,ROUND_MODE);
    }
  else
    {
      mpfr_set_ui(mpfr_A,0,ROUND_MODE);
    }
  at_node->SetMpfr(&mpfr_A,ROUND_MODE);
  for (tav_j=1;tav_j<=n;tav_j++)
    {
      a_new_node=new ListNode(n,tav_j,MY_PREC);
      at_node->SetNext(a_new_node);
      if (!(t_ac_nom==0))
	{
	  function_G(mpfr_A,n,tav_j,t_ac_nom,t_ac_den,ROUND_MODE,MY_PREC);
	}
      else
	{
	  if (tav_j==n)
	    {
	      mpfr_set_ui(mpfr_A,1,ROUND_MODE);
	    }
	  else
	    {
	      mpfr_set_ui(mpfr_A,0,ROUND_MODE);
	    }
	}
      a_new_node->SetMpfr(&mpfr_A,ROUND_MODE);
      at_node=a_new_node;
    }
  a_new_node=NULL;
  at_node=NULL;
  mpfr_clear(mpfr_A);
  mpfr_free_cache();
  return Mpfr_tav_vect;
}

