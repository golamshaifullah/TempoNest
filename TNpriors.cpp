#include <stdio.h>
#include <stdlib.h>
#include "tempo2.h"
#include "TempoNest.h"

//void *globalcontext;

extern void *globalcontext;

void TNprior (double cube[], double theta[], int nDims, void *context)
{
  for(int p=0;p<nDims;p++){
    theta[p]=(((MNStruct *)globalcontext)->PriorsArray[p+nDims]-((MNStruct *)globalcontext)->PriorsArray[p])*cube[p]+((MNStruct *)globalcontext)->PriorsArray[p];
  }

}
