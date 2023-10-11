extern double miyamotonagai_value(double t, double *pars, double *q, int n_dim);
extern void miyamotonagai_gradient(double t, double *pars, double *q, int n_dim, double *grad);
extern void miyamotonagai_hessian(double t, double *pars, double *q, int n_dim, double *hess);
extern double miyamotonagai_density(double t, double *pars, double *q, int n_dim);