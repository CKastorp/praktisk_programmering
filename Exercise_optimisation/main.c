#include"stdio.h"
#include"gsl/gsl_errno.h"
int findmin(void);
int curvefit(void);

int main(){
int status=findmin();
fprintf(stderr,"Findmin status: %i.\n",status);
status=curvefit();
fprintf(stderr,"Curvefit status: %i.\n",status);

return 0;
}
