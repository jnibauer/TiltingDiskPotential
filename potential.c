#include <math.h>


double nan_density(double t, double *pars, double *q, int n_dim) { return NAN; }
double nan_value(double t, double *pars, double *q, int n_dim) { return NAN; }
void nan_gradient(double t, double *pars, double *q, int n_dim, double *grad) {}
void nan_hessian(double t, double *pars, double *q, int n_dim, double *hess) {}

double null_density(double t, double *pars, double *q, int n_dim) { return 0; }
double null_value(double t, double *pars, double *q, int n_dim) { return 0; }
void null_gradient(double t, double *pars, double *q, int n_dim, double *grad){}
void null_hessian(double t, double *pars, double *q, int n_dim, double *hess) {}

/* ---------------------------------------------------------------------------
    Miyamoto-Nagai flattened potential
*/
double miyamotonagai_value(double t, double *pars, double *q, int n_dim) {
    /*  pars:
            - G (Gravitational constant)
            - m (mass scale)
            - a (length scale 1) TODO
            - b (length scale 2) TODO
            - period (time) Yeah thats confusing.
            - n1
            - n2
            - n3
            Following http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf Eq 20
            - Note a minus sign is added. We want the disk to rotate by theta deg around n_hat (counter clockwise). We do this by rotating 
             coords in disk principal axis frame, which involves rotating by the inverse (so -theta deg counterclockwise!)
    */
    double two_pi, thetat, costheta, sintheta, n1, n2, n3, xrot, yrot, zrot, zd;
    two_pi = 2*3.141592653589793238;
    thetat = -(two_pi/pars[4])*t;
    costheta = cos(thetat);
    sintheta = sin(thetat);
    n1 = pars[5];
    n2 = pars[6];
    n3 = pars[7];
    xrot = q[0]*(costheta + (n1*n1)*(1.-costheta)) + q[1]*(n1*n2*(1.-costheta)-n3*sintheta) + q[2]*(n1*n3*(1.-costheta)+n2*sintheta);
    yrot = q[0]*(n1*n2*(1.-costheta) + n3*sintheta) + q[1]*(costheta + (n2*n2)*(1.-costheta)) + q[2]*(n2*n3*(1.-costheta)-n1*sintheta);
    zrot = q[0]*(n1*n3*(1.-costheta)-n2*sintheta) + q[1]*(n2*n3*(1.-costheta)+n1*sintheta) + q[2]*(costheta + (n3*n3)*(1.-costheta));
    
    zd = (pars[2] + sqrt(zrot*zrot + pars[3]*pars[3]));
    return -pars[0] * pars[1] / sqrt(xrot*xrot + yrot*yrot + zd*zd);
}

void miyamotonagai_gradient(double t, double *pars, double *q, int n_dim, double *grad) {
    /*  pars:
            - G (Gravitational constant)
            - m (mass scale)
            - a (length scale 1) TODO
            - b (length scale 2) TODO
            - period (time) Yeah thats confusing.
            - n1
            - n2
            - n3
             Following http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf Eq 20
    */
    double two_pi, thetat, costheta, sintheta, n1, n2, n3, xrot, yrot, zrot, sqrtz, zd, fac, grad0, grad1, grad2;
    two_pi = 2*3.141592653589793238;
    thetat = -(two_pi/pars[4])*t;
    costheta = cos(thetat);
    sintheta = sin(thetat);
    n1 = pars[5];
    n2 = pars[6];
    n3 = pars[7];
    xrot = q[0]*(costheta + (n1*n1)*(1.-costheta)) + q[1]*(n1*n2*(1.-costheta)-n3*sintheta) + q[2]*(n1*n3*(1.-costheta)+n2*sintheta);
    yrot = q[0]*(n1*n2*(1.-costheta) + n3*sintheta) + q[1]*(costheta + (n2*n2)*(1.-costheta)) + q[2]*(n2*n3*(1.-costheta)-n1*sintheta);
    zrot = q[0]*(n1*n3*(1.-costheta)-n2*sintheta) + q[1]*(n2*n3*(1.-costheta)+n1*sintheta) + q[2]*(costheta + (n3*n3)*(1.-costheta));
    

    sqrtz = sqrt(zrot*zrot + pars[3]*pars[3]);
    zd = pars[2] + sqrtz;
    fac = pars[0]*pars[1] * pow(xrot*xrot + yrot*yrot + zd*zd, -1.5);

    grad0 = grad[0] + fac*xrot;
    grad1 = grad[1] + fac*yrot;
    grad2 = grad[2] + fac*zrot * (1. + pars[2] / sqrtz);
    
    
    grad[0] = grad0*(costheta + (n1*n1)*(1.-costheta)) + grad1*(n1*n2*(1.-costheta)+n3*sintheta) + grad2*(n1*n3*(1.-costheta)-n2*sintheta);
    grad[1] = grad0*(n1*n2*(1.-costheta) - n3*sintheta) + grad1*(costheta + (n2*n2)*(1.-costheta)) + grad2*(n2*n3*(1.-costheta)+n1*sintheta);
    grad[2] = grad0*(n1*n3*(1.-costheta)+n2*sintheta) + grad1*(n2*n3*(1.-costheta)-n1*sintheta) + grad2*(costheta + (n3*n3)*(1.-costheta));
        
}

double miyamotonagai_density(double t, double *pars, double *q, int n_dim) {
    /*  pars:
            - G (Gravitational constant)
            - m (mass scale)
            - a (length scale 1) TODO
            - b (length scale 2) TODO
    */

    double M, a, b;
    M = pars[1];
    a = pars[2];
    b = pars[3];

    double R2 = q[0]*q[0] + q[1]*q[1];
    double sqrt_zb = sqrt(q[2]*q[2] + b*b);
    double numer = (b*b*M / (4*M_PI)) * (a*R2 + (a + 3*sqrt_zb)*(a + sqrt_zb)*(a + sqrt_zb));
    double denom = pow(R2 + (a + sqrt_zb)*(a + sqrt_zb), 2.5) * sqrt_zb*sqrt_zb*sqrt_zb;

    return numer/denom;
}

void miyamotonagai_hessian(double t, double *pars, double *q, int n_dim,
                           double *hess) {
    /*  pars:
            - G (Gravitational constant)
            - m (mass scale)
            - a (length scale 1) TODO
            - b (length scale 2) TODO
    */
    double G = pars[0];
    double m = pars[1];
    double a = pars[2];
    double b = pars[3];
    double x = q[0];
    double y = q[1];
    double z = q[2];

    double tmp_0 = pow(x, 2);
    double tmp_1 = pow(y, 2);
    double tmp_2 = pow(z, 2);
    double tmp_3 = pow(b, 2) + tmp_2;
    double tmp_4 = sqrt(tmp_3);
    double tmp_5 = a + tmp_4;
    double tmp_6 = pow(tmp_5, 2);
    double tmp_7 = tmp_0 + tmp_1 + tmp_6;
    double tmp_8 = G*m;
    double tmp_9 = tmp_8/pow(tmp_7, 3.0/2.0);
    double tmp_10 = 3*tmp_8/pow(tmp_7, 5.0/2.0);
    double tmp_11 = tmp_10*x;
    double tmp_12 = -tmp_11*y;
    double tmp_13 = tmp_5/tmp_4;
    double tmp_14 = tmp_13*z;
    double tmp_15 = -tmp_11*tmp_14;
    double tmp_16 = -tmp_10*tmp_14*y;
    double tmp_17 = 1.0/tmp_3;
    double tmp_18 = tmp_2*tmp_9;

    hess[0] = hess[0] + -tmp_0*tmp_10 + tmp_9;
    hess[1] = hess[1] + tmp_12;
    hess[2] = hess[2] + tmp_15;
    hess[3] = hess[3] + tmp_12;
    hess[4] = hess[4] + -tmp_1*tmp_10 + tmp_9;
    hess[5] = hess[5] + tmp_16;
    hess[6] = hess[6] + tmp_15;
    hess[7] = hess[7] + tmp_16;
    hess[8] = hess[8] + -tmp_10*tmp_17*tmp_2*tmp_6 + tmp_13*tmp_9 + tmp_17*tmp_18 - tmp_18*tmp_5/pow(tmp_3, 3.0/2.0);
}
