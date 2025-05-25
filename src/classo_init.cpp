#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Complex.h>

extern "C" {
    void F77_NAME(classocd_warm_)(Rcomplex* x, Rcomplex* y, int* n, int* p, double* lambda, Rcomplex* b0, Rcomplex* b);
    void F77_NAME(classocd_warm_screen_)(Rcomplex* x, Rcomplex* y, int* n, int* p, double* lambda, double* lambda0, Rcomplex* b0, Rcomplex* b);
    void F77_NAME(cglassocd_noscale_)(Rcomplex* s, int* p, double* lambda, Rcomplex* theta, Rcomplex* w, Rcomplex* w0, int* w_init, int* maxiter, double* tol, int* h, int* final_cycle);
    void F77_NAME(cglassocd_scaled_)(Rcomplex* s, int* p, double* lambda, Rcomplex* theta, Rcomplex* w, Rcomplex* w0, int* w_init, int* maxiter, double* tol, int* h, int* final_cycle);
}

static const R_FortranMethodDef FortranEntries[] = {
    // name      pointer      Num args
    {"classocd_warm_", (DL_FUNC)&F77_NAME(classocd_warm_),        7},
    {"classocd_warm_screen_", (DL_FUNC)&F77_NAME(classocd_warm_screen_),        8},
    {"classocd_warm_screen_", (DL_FUNC)&F77_NAME(cglassocd_noscale_),        11},
    {"classocd_warm_screen_", (DL_FUNC)&F77_NAME(cglassocd_scaled_),        11},
    {NULL   ,             NULL,        0}   // Placeholder(?) to indicate last one.
};

void R_init_cxreg(DllInfo* info) {
    R_registerRoutines(
        info,       // DllInfo
        NULL,      // .C
        NULL,      // .Call
        FortranEntries,      // Fortran
        NULL       // External
    );
    R_useDynamicSymbols(info, TRUE);
}