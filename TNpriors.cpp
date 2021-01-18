#include <stdio.h>
#include <stdlib.h>
//#include "TempoNest.h"

//void *globalcontext;
void TNprior (double cube[], double theta[], int nDims, void *context)
{

  for (int i=0; i< nDims; i++) theta[i] = cube[i];
  //for (int i=0; i< nDims; i++) printf("theta[%d] = %lf\n", i, theta[i]);
}
