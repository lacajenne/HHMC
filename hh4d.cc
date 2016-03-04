
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

const int DIM   = 4;  // number of lattice dimensions     
const int DIR_T = 0;  // indices for directions
const int DIR_X = 1;  
const int DIR_Y = 2;
const int DIR_Z = 3;

const int NEW_RUN = 0;  // permits to restart from an old field configuration
const int RESTART = 1;  
const int REJECT  = 0;  // Metropolis 
const int ACCEPT  = 1;

// tables of nearest neighbour sites 
// for all directions(sites have a single index: 0..VOL-1), 
// pos.(p) and neg.(n)
int **nnp; 
int **nnn;

// parameters of the model
// f0 is the initial field value (cold start)
// f : field
// G : phi^2 term
// H : phi^4 term
double G, H, J, f0;         

// lattice extent variables 
int L, VOL;

// unit of momenta (2Pi/L) 
double p_unit;

// random number stuff
long   int idum;
void   init_rg(long int); // random number initialization
double uniform_rg();      // uniform random number in [0,1]
double gaussian_rg();     // gaussian random numbers with zero ave. and unit variance

// HMC functions
double action( double* field );
void   eval_force( double* f_array, double* field );

// lattice geometry
int  coord_to_index( int t, int x, int y, int z );  // find site index from coordinates
void init_lattice(); // assign the table of nearest neighbours

// calculation of some observables
double ave_field(double *field);

// IO utility functions
void SaveLattice(double *field, char *filename);
void LoadLattice(double *field, char *filename);
void Log(char *msg, char *filename);


int main( int argc, char *argv[])
{
  int run_status; // restart? 0 for new run, 1 for restart
  
  int nconf;  //  number of field configurations to be generated
  int nequil; //  number of thermalization field configurations (not used for meas.) 
  
  char logfile[128]; // name of file with results for observables
  char latfile[128]; // name of file with lattice field configuration
  char logbuf[512];  // buffer for results to be piped into logfile

  // correlator data
  FILE *corr_file;
  char file_name_corr[128];

  
  // HMC variables
  double E, En, dE, S, Sn;
  double dt;
  int md_steps;
  int nsave = 10;
  int status;

  // take values of parameters from the command line
  if ( argc != 11 ){
    cerr << "usage : <L> <G> <H> <J> <f0> <nconf> <nequil> <dt> <md_steps> <run_status(new:0, restart:1)>\n";
    exit(-1);
  }
  
  L          = atoi(argv[1]);  
  G          = atof(argv[2]);
  H          = atof(argv[3]);
  J          = atof(argv[4]);
  f0         = atof(argv[5]);
  nconf      = atoi(argv[6]);
  nequil     = atoi(argv[7]);
  dt         = atof(argv[8]);
  md_steps   = atoi(argv[9]);
  run_status = atoi(argv[10]);
 
  // evaluate lattice volume
  VOL = L*L*L*L;
  
  // unit of momenta
  p_unit = 2.0 * M_PI / L;

  // display simulation parameters on screen
  cout << "\nRUN PARAMETERS:\n\n";
  cout << "VOL    : " << L << " x " << L <<  " x " << L << " x " << L << "\n";
  cout << "G      : " << G << "\n";
  cout << "H      : " << H << "\n";
  cout << "J      : " << J << "\n";
  cout << "f0     : " << f0 << "\n";
  cout << "nconf  : " << nconf << "\n"; 
  cout << "nequil : " << nequil << "\n"; 
  cout << "md dt  : " << dt << "\n";
  cout << "md NS  : " << md_steps << "\n";

  if(run_status==NEW_RUN){
    cout << "\nNEW RUN ...\n\n";
  } else {
    cout << "\nRESTART ...\n\n";
  }

  // build name of correlator file
  sprintf(file_name_corr, "corr.dat");
  
  // correct for the decorrelation (see later)
  nconf *= nsave;

  // initialize the lattice geometry and the random numbers
  init_lattice();
  init_rg(-time(NULL));
  
  // the fields
  double *field  = new double[VOL];
  double *nfield = new double[VOL];
  double *mom    = new double[VOL];
  double *momt   = new double[VOL];
  double *force  = new double[VOL];
  double *corr   = new double[L];

  // read field configuration from file or start anew
  if(run_status==NEW_RUN){
    // initialize field
    for(int site=0; site<VOL; site++){
      double r0 = 0.05;
      double r = 1.0 + (2.0*uniform_rg()-1.0)*r0;
      field[site] = f0 * r;
    }
  } else {
    // load last saved field configuration
    sprintf(latfile, "latconf.dat");
    LoadLattice(field, latfile);
  }
  
  for(int site=0; site<VOL; site++) nfield[site] = field[site];
  S = action(field);

  int acc_count = 0;
  
  while( acc_count < nequil + nconf ){
    
    acc_count += 1;

    for(int a=0; a<VOL; a++){
      nfield[a] = field[a];
      mom[a]    = gaussian_rg();
    }
    
    E = S;
    for(int a=0;a<VOL;a++) E += 0.5 * mom[a] * mom[a];
    eval_force(force, nfield);
    for(int a=0;a<VOL;a++) momt[a] = mom[a] + 0.5 * dt * force[a];
    
    for(int step=0;step<md_steps;step++){
      for(int a=0;a<VOL;a++) nfield[a] += dt * momt[a];
      eval_force(force, nfield);
      for(int a=0;a<VOL;a++) momt[a]   += dt * force[a];
    }

    for(int a=0;a<VOL;a++) mom[a] = momt[a] - 0.5 * dt * force[a];
    
    Sn = action( nfield );
    En = Sn;
    for(int a=0;a<VOL;a++) En += 0.5 * mom[a] * mom[a];

    dE = En - E;
    status = REJECT;
    if( dE < 0.0 || uniform_rg() <= exp(-dE) ) status = ACCEPT;

    //cerr << status;
    if(acc_count < nequil){
      sprintf(logbuf, "%.16f", 1.0*status);
      sprintf(logfile, "acc.dat");
      Log(logbuf, logfile);
    }

    if(status == ACCEPT){
      S = Sn;
      for(int a=0;a<VOL;a++) field[a] = nfield[a];
    }
    
    if(acc_count > nequil){
	
      if(acc_count%nsave == 0){
	  
	// average field
	sprintf(logbuf, "%.16f", ave_field(field));
	sprintf(logfile, "phiave.dat");
	Log(logbuf, logfile);

	// correlator
	for(int t=0; t<L; t++) {
	  corr[t] = 0.0;
	  
	  for(int x=0; x<L; x++) {
	    for(int y=0; y<L; y++) {
	      for(int z=0; z<L; z++) {	     
		int site = coord_to_index(t, x, y, z);
		corr[t] += field[site] * field[0];
	      }
	    }
	  }
	  corr[t] /= L*L*L;
	}
	

      }

    }
  }
  
  // save the last field configuration on file
  sprintf(latfile, "latconf.dat");
  SaveLattice(field, latfile);
  
  delete[] field;
  delete[] nfield;
  delete[] mom;
  delete[] momt;
  delete[] force;
  delete[] corr;
  
} // end of main


void init_rg(long int id)
{
  idum = id;
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


double uniform_rg()
{
#ifdef GCC_RAND
  return rand()/((double)RAND_MAX);
#else
  
  //tratto da http://www.library.cornell.edu/nr/cbookcpdf.html 
  
  int j;
  long int k;
  static long int idum2 = 123456789;
  static long int iy = 0;
  static long iv[NTAB];
  double temp;
  
  if(idum <= 0){
    if(-idum < 1)
      idum = 1;
    else
      idum = -idum;
    idum2 = idum;
    for(j = NTAB+7; j >= 0; j--){
      k = idum/IQ1;
      idum = IA1*(idum-k*IQ1)-k*IR1;
      if(idum < 0) idum += IM1;
      if(j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
  }
  k = idum/IQ1;
  idum = IA1*(idum-k*IQ1)-k*IR1;
  if(idum < 0) idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if(iy < 1) iy += IMM1;
  if((temp=AM*iy) > RNMX) return RNMX;
  return temp;
#endif
}

double gaussian_rg()
{
  // Genera x con distribuzione exp(-.5*x^2)dx in (-\infty,+\infty).                                                                                     
  
  double v1, v2, zsq, R;
  static int flag = 0;
  static double gset;
  if(flag == 0){
    flag = 1;
    do{
      v1 = 2.0*uniform_rg()-1.0;
      v2 = 2.0*uniform_rg()-1.0;
      zsq = v1*v1 + v2*v2;
    } while(zsq > 1 || zsq == 0);
    R = sqrt(-2.*log(zsq)/zsq);
    gset = v1 * R;
    return v2 * R;
  }

  flag = 0;
  return gset;
}


inline int coord_to_index( int t, int x, int y, int z )
{
  // return the site index from the coordinates (t, x, y) 
  return t + L*x + L*L*y + L*L*L*z;
}


void init_lattice() {
  
  nnp = new int*[VOL];
  nnn = new int*[VOL];
  
  for(int a=0;a<VOL;a++){
    
    nnp[a] = new int[DIM]; // one n.n for each pos. dir. for each lattice point (a)
    nnn[a] = new int[DIM]; // one n.n for each neg. dir. for each lattice point (a)
    
    // extract the coordinates (t, x, y, z) from the site index
    int t, x, y, z, b;
    t = a % L;
    b = a / L;
    x = b % L;
    b = a / (L*L);
    y = b % L;
    b = a / (L*L*L);
    z = b % L;
    
    // identify nearest neighbours implementing periodic b.c.

    nnp[a][DIR_T] = t == (L-1) ? coord_to_index(0, x, y, z)   : coord_to_index(t+1, x, y, z);
    nnn[a][DIR_T] = t == 0     ? coord_to_index(L-1, x, y, z) : coord_to_index(t-1, x, y, z);

    nnp[a][DIR_X] = x == (L-1) ? coord_to_index(t, 0, y, z)   : coord_to_index(t, x+1, y, z);
    nnn[a][DIR_X] = x == 0     ? coord_to_index(t, L-1, y, z) : coord_to_index(t, x-1, y, z);

    nnp[a][DIR_Y] = y == (L-1) ? coord_to_index(t, x, 0, z)   : coord_to_index(t, x, y+1, z);
    nnn[a][DIR_Y] = y == 0     ? coord_to_index(t, x, L-1, z) : coord_to_index(t, x, y-1, z);

    nnp[a][DIR_Z] = z == (L-1) ? coord_to_index(t, x, y, 0)   : coord_to_index(t, x, y, z+1);
    nnn[a][DIR_Z] = z == 0     ? coord_to_index(t, x, y, L-1) : coord_to_index(t, x, y, z-1);

  }
}

double action( double* field )
{
  double S = 0.0;
  double f;
  double tb;

  for(int a=0; a<VOL; a++){
    
    f = field[a];
    
    S += DIM * f * f;
    for(int mu=0; mu<DIM; mu++) S -= f * field[nnp[a][mu]];
    
    S += 0.5*G*f*f + 0.25*H*f*f*f*f + J*f;
    
  }

  return S;
}
    

void eval_force( double* f_array, double* field )
{
  double f;

  for(int a=0; a<VOL; a++){
    f = field[a];
    f_array[a] = -2*DIM*f - G*f - H*f*f*f - J;
    for(int mu=0;mu<DIM;mu++){
      f_array[a] += field[nnp[a][mu]] + field[nnn[a][mu]];
    }   
  }
}


inline double ave_field(double *field)
{
  double ret = 0.0;
  for(int site=0; site<VOL; site++) ret += field[site];
  ret /= VOL;
  return ret;
}

inline void SaveLattice(double *field, char *filename)
{
  FILE *file;
  file = fopen(filename, "w");
  for(int i=0;i<VOL;i++){
    fprintf(file, "%.40f\n", field[i]);
  }
  fflush(file);
  fclose(file);
}


inline void LoadLattice(double *field, char *filename)
{
  FILE *file;
  file = fopen(filename, "r");
  for(int i=0;i<VOL;i++){
    int p = 0;
    p = fscanf(file, "%lf\n", &field[i]);
  }
  fflush(file);
  fclose(file);
}


inline void Log(char *msg, char *filename){
  FILE *file;
  file = fopen(filename, "a");
  fprintf(file, "%s\n", msg);
  fflush(file);
  fclose(file);
}

