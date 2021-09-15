
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
using namespace std;
#include "header_files/main_prob_funcs_new.h"
#include "header_files/linked_mpfr_lists.h"
#include "header_files/generate_Tav_list.h"


int main(int argc, char *argv[])
{

  ListNode* Mpfr_tav_vect=NULL;
  ListNode* temp_node=NULL;
  mpfr_t* an_mpfr_ptr=NULL;

  mpfr_t mpfr_prob;
  mp_rnd_t ROUND_MODE=GMP_RNDN; //GMP_RNDU;GMP_RNDD;GMP_RNDN
  mp_prec_t MY_PREC;

  int i;
  unsigned long int n,t_ac_nom,t_ac_den,temp_ui_i;
  long double the_out_val;

  int arg_i, temp_int;

  vector<string> args;
  for (arg_i=1; arg_i<argc; arg_i++)
    {
      args.push_back(argv[arg_i]);
    }
  if (args.size()==4)
    {
      temp_int=atoi(args[0].c_str());
      n=static_cast<unsigned long int>(temp_int);
      temp_int=atoi(args[1].c_str());
      t_ac_nom=static_cast<unsigned long int>(temp_int);
      temp_int=atoi(args[2].c_str());
      t_ac_den=static_cast<unsigned long int>(temp_int);
      MY_PREC=atoi(args[3].c_str());
    }
  else
    {
      cout<<"number of arguments is "<<args.size()<<"\n";
      cout<<"Number of arguments should be 4!!\n";
    }

  Mpfr_tav_vect=function_Tav_1D(n,t_ac_nom,t_ac_den, ROUND_MODE,MY_PREC);
  mpfr_init2(mpfr_prob,MY_PREC);
  cout.precision(20);
  for (i=0;i<n;i++)
    {
      temp_ui_i=n-i;
      temp_node=Mpfr_tav_vect->GetBetaAlphaNode(Mpfr_tav_vect,n,temp_ui_i);
      an_mpfr_ptr=temp_node->GetMpfr();
      mpfr_set(mpfr_prob,*an_mpfr_ptr,ROUND_MODE);
      temp_node=NULL;
      an_mpfr_ptr=NULL;
      the_out_val=mpfr_get_ld(mpfr_prob,ROUND_MODE);
      cout<<temp_ui_i<<"\t"<<the_out_val<<endl;
    }
  mpfr_clear(mpfr_prob);
  Mpfr_tav_vect->EraseList(Mpfr_tav_vect);
  mpfr_free_cache();
  return 0;
}

