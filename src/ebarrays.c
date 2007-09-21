#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>



SEXP makeProbVect(SEXP zArg, SEXP pArg)
{
     double *z, *ansptr;
     double *p;
     SEXP dim;
     SEXP ans;
     int m, n, i, j, k;

     z = NUMERIC_DATA(zArg);
     p = NUMERIC_DATA(pArg);
     dim = GET_DIM(zArg);
     m = INTEGER(dim)[0];
     n = INTEGER(dim)[1];

     /*printf("\nDimensions of z: %d  %d", m, n);*/
     PROTECT(ans = NEW_NUMERIC(m * n));
     ansptr = NUMERIC_POINTER(ans);
     SET_DIM(ans, dim);

     for (i = 0; i < m; i++)
         for (j = 0; j < n; j++)
	 {
	     ansptr[m * j + i] = 0;
	     for (k = 0; k < n; k++)
		 ansptr[m * j + i] += exp(z[m * k + i] + log(p[k]) - z[m * j + i] - log(p[j]));
	     ansptr[m * j + i] = 1 / ansptr[m * j + i];
	 }

     UNPROTECT(1);
     return ans;
}

SEXP makeProbVectArr(SEXP zArg, SEXP pArg)
{
     double *z, *ansptr;
     double *p;
     SEXP dim;
     SEXP ans;
     int m, n, i, j, k, l;

     z = NUMERIC_DATA(zArg);
     p = NUMERIC_DATA(pArg);
     dim = GET_DIM(zArg);
     m = INTEGER(dim)[0];
     n = INTEGER(dim)[1];
     l = INTEGER(dim)[2];

     /*printf("\nDimensions of z: %d  %d %d", m, n, p);*/
     PROTECT(ans = NEW_NUMERIC(m * n * l));
     ansptr = NUMERIC_POINTER(ans);
     SET_DIM(ans, dim);

     for (i = 0; i < m; i++)
         for (j = 0; j < (n*l); j++)
	 {
	     ansptr[m * j + i] = 0;
	     for (k = 0; k < (n*l); k++)
		 ansptr[m * j + i] += exp(z[m * k + i] + log(p[k]) - z[m * j + i] - log(p[j]));
	     ansptr[m * j + i] = 1 / ansptr[m * j + i];
	 }

     UNPROTECT(1);
     return ans;
}
