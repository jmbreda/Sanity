#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <omp.h>
#include <time.h>

using namespace std;

// Compile :
// g++ -std=c++11 -O2 -ffast-math -O3 -fopenmp compute_gene_correlations.cpp -o Sanity_gene_correlations

// Distances functions :
// with error bars:

double get_correlation_simple(double *deli,double *delj,int N_cell);
void print_results(string out_file, double **Rs, double **Rfull, double **Rerr, int *orignum, int N_gene);
double deriv(double r, double *cov,double *phi, double *V, int N_cell);
void get_correlation_full(double *rfull, double *rerr, double *deli, double *eps2i, double vi, double *delj, double *eps2j, double vj, double rs, int N_cell);

// I/O
void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, int &N_gene, int &N_cell);
static void show_usage(void);
void load_matrix(string sanity_folder, double **delta, double **epsilon2, double *variance, int N_gene, int N_cell);


int main(int argc, char** argv){  
	// Declare arguments :
	string sanity_folder("");
	double s2n_cutoff = 1.0;
	int N_gene(0);
	int N_cell(0);
	parse_argv(argc, argv, sanity_folder, s2n_cutoff, N_gene, N_cell);

	cout << "Sanity folder: " << sanity_folder << "\n";
	cout << "Signal/noise cut-off: " << s2n_cutoff << "\n";
	cout << "Nr. of gene: " << N_gene << "\n";
	cout << "Nr. of cell: " << N_cell << "\n";

	// Declare temporary variables
	int g, c;
	double **delta_all = new double *[N_gene];
	for(g=0;g<N_gene;g++)
	  {
	    delta_all[g] = new double[N_cell];
	  }

	double **epsilon2_all = new double *[N_gene];
	for(g=0;g<N_gene;g++)
	  epsilon2_all[g] = new double[N_cell];
	
	double *variance_all = new double [N_gene];

	// Get data
	cout << "Getting data... ";
	load_matrix(sanity_folder, delta_all, epsilon2_all, variance_all, N_gene, N_cell);
	cout << "data loaded \n";
      

	// Compute signal 2 noise
	double s2n_gene[N_gene];
	int N_gene_cutoff(0);
	if (s2n_cutoff > 0.0)
	  {
	    double mean_delta;
	    double var_delta;
	    double mean_epsilon2;
	    
	    for(g=0;g<N_gene;g++)
	      {
		s2n_gene[g] = 0;
		for(c=0;c<N_cell;c++)
		  {
		    s2n_gene[g] += delta_all[g][c]*delta_all[g][c]/epsilon2_all[g][c];
		  }
	        s2n_gene[g] /= ((double) N_cell);
		s2n_gene[g] = sqrt(s2n_gene[g]);
		
		if(s2n_gene[g] >= s2n_cutoff)
		  {
		    N_gene_cutoff += 1;
		  }
	      }
	  }
	else
	  {
	    N_gene_cutoff = N_gene;
	  }

	cout << "Done s2n determination. Number of genes left " << N_gene_cutoff << "\n";

	
	double **delta = new double *[N_gene_cutoff];
	for(g=0;g<N_gene_cutoff;g++)
	  delta[g] = new double[N_cell];

	double **epsilon2 = new double *[N_gene_cutoff];
	for(g=0;g<N_gene_cutoff;g++)
	  epsilon2[g] = new double[N_cell];

	double *variance = new double [N_gene_cutoff];
	int *orignum = new int [N_gene_cutoff];

	if (s2n_cutoff > 0.0)
	  {
	    // Copy genes with signal to noise above cut-off
	    int k = 0;
	    for(g=0;g<N_gene;g++)
	      {
		if(s2n_gene[g] >= s2n_cutoff)
		  {
		    for(c=0;c<N_cell;c++)
		      {
			delta[k][c] = delta_all[g][c];
			epsilon2[k][c] = epsilon2_all[g][c];
		      }
		    variance[k] = variance_all[g];
		    orignum[k] = g;
		    k++;
		  }
	      }
	}
	else
	  {
	    for(g=0;g<N_gene;g++)
	      {
		for(c=0;c<N_cell;c++)
		  {
		    delta[g][c] = delta_all[g][c];
		    epsilon2[g][c] = epsilon2_all[g][c];
		  }
		variance[g] = variance_all[g];
		orignum[g] = g;
	      }
	  }
	
	// delete temporary variables
	for(g=0;g<N_gene;g++)
	  {
	    delete[] delta_all[g];
	    delete[] epsilon2_all[g];
	  }
	delete[] delta_all;
	delete[] epsilon2_all;
	delete[] variance_all;
	

	N_gene = N_gene_cutoff;
	//Initialize space for simple and full correlations
	cout << "Calculating correlations... ";
	double **Rs = new double *[N_gene];
	double **Rfull = new double *[N_gene];
	double **Rerr = new double *[N_gene];
	for(g=0; g<N_gene; g++)
	  {
	    Rs[g] = new double[N_gene];
	    Rfull[g] = new double[N_gene];
	    Rerr[g] = new double [N_gene];
	  }


	//loop over all pairs to get simple correlations first
	int i,j;
	for(i=0;i<N_gene;i++)
	  {
	    for(j=(i+1);j<N_gene;j++)
	      {
		Rs[i][j] = get_correlation_simple(delta[i],delta[j],N_cell);
		//cout << " genes " << i << " and " << j << " Rsimple " << Rs[i][j] << "\n";
	      }
	  }
	cout << "done with simple correlations \nDoing full correlations with errobars...";

       
	//Now first rescale all the deltas and epsilons
	double factor;
	for(g=0;g<N_gene;g++)
	  {
	    for(c=0;c<N_cell;c++)
	      {
		if(variance[g]-epsilon2[g][c] > 0.0)
		  factor = variance[g]/(variance[g]-epsilon2[g][c]);
		else
		  factor = variance[g]/0.000001;
		
		delta[g][c] = delta[g][c]*factor;
		epsilon2[g][c] = epsilon2[g][c]*factor;
	      }
	  }

	//Now do the full correlation calculations
	for(i=0;i<N_gene;i++)
	  {
	    //cout << "gene number " << i << "\n";
	    for(j=(i+1);j<N_gene;j++)
	      {
		get_correlation_full(&(Rfull[i][j]),&(Rerr[i][j]),delta[i],epsilon2[i],variance[i],delta[j],epsilon2[j],variance[j],Rs[i][j],N_cell);
	      }
	  }
	cout << "done\n";
	
	// Print correlations
	string out_file = sanity_folder + "gene_correlations";
	
	if(s2n_cutoff > 0.0)
	  {
	    string s2n_str = to_string(s2n_cutoff);
	    while(s2n_str.back() == '0' || s2n_str.back()=='.')
	      s2n_str = s2n_str.substr(0, s2n_str.size()-1);
	    out_file += "_s2n_gt_" + s2n_str;
	  }
	out_file += ".txt";

	cout << "Printing results in " << out_file;
	print_results(out_file,Rs,Rfull,Rerr,orignum,N_gene);
  	cout << "\n";
	
	return 0;
}

//functions

double get_correlation_simple(double *deli,double *delj, int N_cell){
  int c;
  double mui, muj, vi, vj, cov, r;

  mui = 0;
  muj = 0;
  vi = 0;
  vj = 0;
  cov = 0;
  for(c=0;c<N_cell;++c)
    {
      mui += deli[c];
      muj += delj[c];
      vi += deli[c]*deli[c];
      vj += delj[c]*delj[c];
      cov += deli[c]*delj[c];
    }
  mui /= ((double) N_cell);
  muj /= ((double) N_cell);
  vi /= ((double) N_cell);
  vj /= ((double) N_cell);
  cov /= ((double) N_cell);
  vi -= mui*mui;
  vj -= muj*muj;
  cov -= mui*muj;
  r = cov/sqrt(vi*vj);

  return r;
}


double deriv(double r, double *cov,double *phi, double *V, int N_cell){
  int c;
  double deriv = 0.0;
  double denom;

  for(c=0;c<N_cell;++c)
    {
      denom = 1.0/(phi[c]-r*r);
      deriv += r*denom;
      deriv += cov[c]*denom;
      deriv -= r*(V[c]-2.0*r*cov[c])*denom*denom;
    }

  return deriv;
}


void get_correlation_full(double *rfull, double *rerr, double *deli, double *eps2i, double vi, double *delj, double *eps2j, double vj, double rs, int N_cell){

  //These arrays hold simplified stats for the two genes
  static double *phi = new double[N_cell];
  static double *V = new double[N_cell];
  static double *cov = new double[N_cell];

  int c;
  double sigi = sqrt(vi);
  double sigj = sqrt(vj);
  //Calculate the arrays
  for(c=0;c<N_cell;++c)
    {
      double phi_i = 1.0 + eps2i[c]/vi;
      double phi_j = 1.0 + eps2j[c]/vj;
      double zi = deli[c]/sigi;
      double zj = delj[c]/sigj;

      cov[c] = zi * zj;
      phi[c] = phi_i * phi_j;
      V[c] = zi*zi*phi_j + zj*zj*phi_i;
    }
  
  //Get the derivative of the log-likelihood at the simple guess rs
  double ders = deriv(rs,cov,phi,V,N_cell);

  double rmin, rmax,der;
  if(ders > 0)
    {
      rmin = rs;
      rmax = 1.0;

      der = deriv(rmax,cov,phi,V,N_cell);
      if(der > 0)
	{
	  //Optimum is at r=1. We estimate error-bar by by slope of log-lik at r=1.0 (i.e. assume exponential dist around r=1)
	  *rfull = 1.0;
	  *rerr = -log(1.0-0.159)/der;//note 0.156 is z-score at 1 sigma of a Gaussian. This is fairly arbitrary
	  return;
	}
    }
  else
    {
      rmax = rs;
      rmin = -1.0;

      der = deriv(rmin,cov,phi,V,N_cell);
      if(der < 0)
	{
	  *rfull = -1.0;
	  *rerr = log(1-0.159)/der;
	  return;
	}
    }

  double rtol = 0.0025;
  double curder, rmid;
  while(rmax-rmin > rtol)
    {
      //derivative at the middle
      rmid = 0.5*(rmin+rmax);
      curder = deriv(rmid,cov,phi,V,N_cell);

      if(curder > 0)
	rmin = rmid;
      else
	rmax = rmid;
    }

  double rfinal = 0.5*(rmin+rmax);
  *rfull = rfinal;

  //Now calculate the error bar (expression is a bit messy)
  double secder = 0;
  double num, denom;
  for(c=0;c<N_cell;++c)
    {
      num = phi[c]*phi[c] + 6.0 *cov[c]*phi[c]*rfinal+2.0*cov[c]*rfinal*rfinal*rfinal-rfinal*rfinal*rfinal*rfinal-V[c]*(phi[c]+3.0*rfinal*rfinal);
      denom = 1.0/(phi[c]-rfinal*rfinal);
      secder += num*denom*denom*denom;
    }
  *rerr = 1/sqrt(-secder);

  return;
}

void load_matrix(string sanity_folder, double **delta, double **epsilon2, double *variance, int N_gene, int N_cell){

	string my_file;
	FILE *infp;
	char ss[3000000];
	char *token;
	char *retval;
	int g, c;
	double epsilon;

	// Read delta
	my_file = sanity_folder + "delta.txt";
	infp = (FILE *) fopen(my_file.c_str(),"r");
	if(infp == NULL)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",my_file.c_str());
	    exit(EXIT_FAILURE);
	  }
	for (g = 0; g < N_gene; g++)
	  {
	    retval = fgets(ss,3000000,infp);
	    if(retval == NULL)
	      {
		fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
		return;
	      }
	    token = strtok(ss," \t");
	    for (c=0; c<N_cell; c++)
	      {	
		if(token != NULL)
		  {
		    delta[g][c] = atof(token);
		  }
		else
		  {
		    fprintf(stderr,"Error: not enough fields on line number %d:\n",g);
		    return; 
		  }
		token = strtok(NULL," \t");
	      }
	    token = strtok(NULL," \t");
	    if(token != NULL)
	      {
		fprintf(stderr,"Error: too many fields on line number %d:\n",g);
	      }
	  }
	fclose(infp);

	// Read epsilon
	my_file = sanity_folder + "d_delta.txt";
	infp = (FILE *) fopen(my_file.c_str(),"r");
	if(infp == NULL)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",my_file.c_str());
	    exit(EXIT_FAILURE);
	  }
	for (g = 0; g < N_gene; g++)
	  {
	    retval = fgets(ss,3000000,infp);
	    if(retval == NULL)
	      {
		fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
		return;
	      }
	    token = strtok(ss," \t");
	    for (c=0; c<N_cell; c++)
	      {	
		if(token != NULL){
		  epsilon = atof(token);
		  epsilon2[g][c] = epsilon*epsilon;
		}
		else
		  {
		    fprintf(stderr,"Error: not enough fields on line number %d:\n",g);
		    return; 
		  }
		token = strtok(NULL," \t");
	      }
	    token = strtok(NULL," \t");
	    if(token != NULL)
	      {
		fprintf(stderr,"Error: too many fields on line number %d:\n",g);
	      }
	  }
	fclose(infp);

	
	my_file = sanity_folder + "variance.txt";
	infp = (FILE *) fopen(my_file.c_str(),"r");
	if(infp == NULL)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",my_file.c_str());
	    exit(EXIT_FAILURE);
	  }
	for (g = 0; g < N_gene; g++)
	  {
	    retval = fgets(ss,3000000,infp);
	    if(retval == NULL)
	      {
		fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
		return;
	      }
	    token = strtok(ss," \t");
	    if(token != NULL)
	      {
		variance[g] = atof(token);
	      }
	    else
	      {
		fprintf(stderr,"Error: not enough fields on line number %d:\n",g);
		return;
	      }
	    token = strtok(NULL," \t");
	    if(token != NULL)
	      {
		fprintf(stderr,"Error: too many fields on line number %d:\n",g);
	      }
	  }
	fclose(infp);

	return;
}

void print_results(string out_file, double **Rs, double **Rfull, double **Rerr, int *orignum, int N_gene){

  	
	int i, j;
	string file_name(out_file);
	ofstream outfile(file_name);
	for (i=0; i<(N_gene-1); i++)
	  {
	    for (j=i+1; j<N_gene; j++)
	      {
		outfile << orignum[i] << "\t" << orignum[j] << "\t" << Rs[i][j] << "\t" << Rfull[i][j] << "\t" << Rerr[i][j] << "\n";
	      }
	  }
	outfile.close();
}


void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, int &N_gene, int &N_cell){

    if (argc<2)
        show_usage();

    int i;

    string get_help [2];
    get_help[0] = "-h";
    get_help[1] = "--help";

    for(i=1;i<argc;i++){
        if (argv[i] == get_help[0] || argv[i] == get_help[1])
            show_usage();
    }

    string get_version [2];
    get_version[0] = "-v";
    get_version[1] = "--version";
    for(i=1;i<argc;i++)
      {
        if (argv[i] == get_version[0] || argv[i] == get_version[1])
	  {
            cout << "v1.0" << "\n";
            exit(0);
	  }
      }

    int N_param(2);
    string error_bar("true");
    string to_find[2][2] = {{"-f", "--folder"},
                            {"-s2n", "--signal_to_noise_cutoff"}};
    
    int j;
    int idx;
    for(j=0;j<N_param;j++){
        idx = 0;
        for(i=1;i<argc;i++){
            if (argv[i] == to_find[j][0] || argv[i] == to_find[j][1]){
                idx = i;
                if ( idx+1 > argc-1 ){
                    cerr << "Error in argument parsing :\n"
                         << argv[i] << " option missing\n";
                    show_usage();
                }
                if(j==0) sanity_folder = argv[idx+1];
                if(j==1) s2n_cutoff = atof(argv[idx+1]);
		
                // add '/' to out_folder if not already
                if( j == 0 && sanity_folder.back() != '/' )
                    sanity_folder = sanity_folder + '/';
            }
        }
        if (idx == 0 && j == 0){
            cerr << "Error in argument parsing :\n"
            << "missing input folder name\n";
            show_usage();
        }
    }
    
    // Get number of gene and number of cells
    N_gene = -1;
    N_cell = -1;
    ifstream in(sanity_folder + "delta.txt");
    if (in.is_open()){
      N_gene = 0;
      N_cell = 0;
      string line;
      
      // get first line and count columns
      getline(in, line);
      ++N_gene;
      istringstream iss(line);
      do{
	string sub;
	iss >> sub;
	if (sub.length())
	  ++N_cell;
      }
      while(iss);
      
      // Now count genes in remaining lines
      while ( getline(in, line) ){
	++N_gene;
      }
    }else{
      cerr << "Unable to open " + sanity_folder + "delta.txt\n";
    }

    return;
    
}

static void show_usage(void)
{
    cerr << "Usage: Sanity_gene_correlations <option(s)>\n"
         << "Options:\n"
         << "\t-h,--help\t\tShow this help message\n"
         << "\t-v,--version\t\tShow the current version\n"
         << "\t-f,--folder\t\tSpecify the input folder with extended output from Sanity\n"
         << "\t-s2n,--signal_to_noise_cutoff\tMinimal signal/noise of genes to include in the distance calculation (default: 1.0)\n";
    exit(0);
}

