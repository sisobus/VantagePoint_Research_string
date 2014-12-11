#define MYPI  3.141592653 
#include "SDPdef.h" 
 
typedef struct { 
     int    index;  /* index of value; needed after sorting the values */ 
     double val;    /* the vector of doubles                           */ 
 } orderCut; 
 
orderCut *OrdAlloc(int  len, char *info); 
void OrdFree(orderCut **x);
void GetSinve(sdpdat *,double *);
void GetXrow( sdpdat *pd,double *sinve,double *xrowi,double *ei );
double EvalObjective( symat, double *, int, int); 
void DefineCut(double  *, int, orderCut *, int,int,int ); 
void GetRandVec(double *, int); 
void RandReduce( sdpdat * ); 
void SaveIt(double *, int, int, char *);
 
