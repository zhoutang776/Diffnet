#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void diffnet_lasso(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void diffnet_mcp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void diffnet_scad(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"diffnet_lasso", (DL_FUNC) &diffnet_lasso, 16},
    {"diffnet_mcp",   (DL_FUNC) &diffnet_mcp,   16},
    {"diffnet_scad",  (DL_FUNC) &diffnet_scad,  16},
    {NULL, NULL, 0}
};

void R_init_diffnet(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}