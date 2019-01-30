typedef struct {
	int systemFlag;
	double size;
	double dr;
	int ensembleFlag;
	int N;
	double M;
	double Ma;
	double den, den_v;
	double pl, pv;

	double mu;
	double R;
	double Ri;
	double RSolute;

	double step;

	double HSDiam;
	double sig;
	double T;
	double T_assoc;

	double rCut;
	
	int FHSFlag;
	int FDispFlag;

	double q;
	double tol;


	int wallFlag;
	double wallT;
} INPUT;
