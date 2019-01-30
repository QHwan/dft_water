typedef struct {
	int systemFlag;
	double size;
	double dr;
	int ensembleFlag;
	int N;
	double den, den_v;

	double mu;
	double Ri;

	double step;

	double HSDiam;
	double sig;
	double T;

	double rCut;
	
	int FHSFlag;
	int FDispFlag;

	double q;
	double tol;

	int pressure;
} INPUT;
