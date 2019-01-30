
double get_a(double T,double rc);

double fmu_id(double rho);

double fp_id(double rho);

double fmu_hs(double rho);

double fp_hs(double rho);

double fmu_disp(double rho, double a);

double fp_disp(double rho, double a);

double get_ghs_bulk(double rho);

double get_gphs_bulk(double rho);

double get_D_bulk(double rho,double K,double T_assoc);

double get_Dp_bulk(double rho,double K,double T_assoc);

double get_X_bulk(double rho,double Ma,double K,double T_assoc);

double get_Xp_bulk(double rho,double Ma,double K,double T_assoc);

double fmu_assoc(double rho,double Ma,double K,double T_assoc);

double fp_assoc(double rho,double Ma,double K,double T_assoc);


void find_roots_mu(double *dVec, double T, double a, double mu_test);


int checkTol(double *tolVec, double *denVec, double *denOutVec, int nGrids, double tol);
double returnTol(double *tolVec, double *denVec, double *denOutVec, int nGrids, double tol);


void updateDen(double *denVec, double *denInVec, double *denOutVec, int nGrids, double q);
