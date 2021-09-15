#ifndef LINKED_LISTS_MPFR
#define LINKED_LISTS_MPFR

class ListNode
{
private:
  unsigned long int MyBeta;
  unsigned long int MyAlpha;
  ListNode* MyNext;
  mpfr_t MyMpfr;
public:
  ListNode(unsigned long int a_beta,unsigned long int an_alpha,mp_prec_t MY_PREC);//Constructor
  ~ListNode();//Destructor
  void SetMpfr(mpfr_t* mpfr_ptr,mp_rnd_t ROUND_MODE){mpfr_set(MyMpfr,*mpfr_ptr,ROUND_MODE);}
  void SetNext(ListNode* a_node){MyNext=a_node;}
  unsigned long int GetBeta(){return MyBeta;}
  unsigned long int GetAlpha(){return MyAlpha;} 
  mpfr_t* GetMpfr(){return &MyMpfr;}
  ListNode* GetNext(){return MyNext;}
  ListNode* GetBetaAlphaNode(ListNode* a_list,unsigned long int a_beta,unsigned long int an_alpha);
  void EraseList(ListNode* a_list);
  void DebugPrintList(ListNode* a_list,mp_rnd_t ROUND_MODE);
};

#endif

ListNode::ListNode(unsigned long int a_beta,unsigned long int an_alpha,mp_prec_t MY_PREC)
{
  MyBeta=a_beta;
  MyAlpha=an_alpha;
  mpfr_init2(MyMpfr,MY_PREC);
  MyNext=NULL;
}
ListNode::~ListNode()
{
  MyBeta=0;
  MyAlpha=0;
  mpfr_clear(MyMpfr);
  MyNext=NULL;
}
ListNode* ListNode::GetBetaAlphaNode(ListNode* a_list,unsigned long int a_beta,unsigned long int an_alpha)
{
  unsigned long int at_beta,at_alpha;
  ListNode* at_node;
  at_node=a_list;
  at_beta=at_node->GetBeta();
  while (at_beta<a_beta)
    {
      at_node=at_node->GetNext();
      at_beta=at_node->GetBeta();
    }
  at_alpha=at_node->GetAlpha();
  while (at_alpha<an_alpha)
    {
      at_node=at_node->GetNext();
      at_alpha=at_node->GetAlpha();
    }
  return at_node;
}
void ListNode::EraseList(ListNode* a_list)
{
  ListNode* next_node;
  while (!(a_list==NULL))
    {
      next_node=a_list->GetNext();
      delete a_list;
      a_list=next_node;
    }
}
void ListNode::DebugPrintList(ListNode* a_list, mp_rnd_t ROUND_MODE)
{
  mpfr_t* an_mpfr_ptr;
  ListNode* next_node;
  long double out_val;
  while (!(a_list==NULL))
    {
      an_mpfr_ptr=a_list->GetMpfr();
      out_val=mpfr_get_ld(*an_mpfr_ptr,ROUND_MODE);
      cout<<"Beta="<<a_list->GetBeta()<<"\tAlpha="<<a_list->GetAlpha()<<"\tvalue="<<out_val<<endl;
      next_node=a_list->GetNext();
      a_list=next_node;
    }
}




