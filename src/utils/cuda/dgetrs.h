/*!**********************************************************************
 * \file dgetrs.h
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-22.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DGETRS_H_ZG149ZBV
#define DGETRS_H_ZG149ZBV

int dgetrscuda_(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
	
int dtrsmcuda_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
	
int xerblacuda_(char *srname, int *info);
	
bool lsamecuda_(char *ca, char *cb);

#endif /* end of include guard: DGETRS_H_ZG149ZBV */
