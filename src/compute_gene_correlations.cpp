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

#include "Version.h"

using namespace std;

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// The SI does not describe the gene-gene correlations or the signal-to-noise cut-off, so only delta, epsilon and v_g below refer to it.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells.
// Two letters mean a matrix stored in that order, so log_fold_change_gc is indexed [gene][cell] and corr_simple_gg is genes x genes.
// A suffix such as _i_c is a vector over the cells for gene i.

// Compile :
// g++ -std=c++11 -O2 -ffast-math -O3 -fopenmp compute_gene_correlations.cpp -o Sanity_gene_correlations

// Distances functions :
// with error bars:

double get_correlation_simple(double *log_fold_change_i_c,double *log_fold_change_j_c,int num_cells);
void print_results(string out_file, double **corr_simple_gg, double **corr_full_gg, double **corr_error_gg, int *original_gene_index_g, int num_genes);
double log_likelihood_slope(double corr, double *z_product_c,double *variance_inflation_product_c, double *cross_term_c, int num_cells);
void get_correlation_full(double *corr_full, double *corr_error, double *log_fold_change_i_c, double *log_fold_change_error_sq_i_c, double gene_variance_i, double *log_fold_change_j_c, double *log_fold_change_error_sq_j_c, double gene_variance_j, double corr_simple, int num_cells);

// I/O
void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, int &num_genes, int &num_cells);
static void show_usage(void);
void load_matrix(string sanity_folder, double **log_fold_change_gc, double **log_fold_change_error_sq_gc, double *gene_variance_g, int num_genes, int num_cells);


int main(int argc, char** argv){  
	// Declare arguments :
	string sanity_folder("");
	double s2n_cutoff = 1.0;
	int num_genes(0);
	int num_cells(0);
	parse_argv(argc, argv, sanity_folder, s2n_cutoff, num_genes, num_cells);

	cout << "Sanity folder: " << sanity_folder << "\n";
	cout << "Signal/noise cut-off: " << s2n_cutoff << "\n";
	cout << "Nr. of gene: " << num_genes << "\n";
	cout << "Nr. of cell: " << num_cells << "\n";

	// Declare temporary variables
	int g, c;
	double **log_fold_change_all_gc = new double *[num_genes]; // <delta_gc> (SI Eq. 39), indexed [gene][cell], read from delta.txt; before the signal-to-noise filter.
	for(g=0;g<num_genes;g++)
	  {
	    log_fold_change_all_gc[g] = new double[num_cells];
	  }

	double **log_fold_change_error_sq_all_gc = new double *[num_genes]; // epsilon_gc^2 (SI Eq. 42), from the error bars in d_delta.txt.
	for(g=0;g<num_genes;g++)
	  log_fold_change_error_sq_all_gc[g] = new double[num_cells];
	
	double *gene_variance_all_g = new double [num_genes]; // v_g (SI Eq. 17), read from variance.txt.

	// Get data
	cout << "Getting data... ";
	load_matrix(sanity_folder, log_fold_change_all_gc, log_fold_change_error_sq_all_gc, gene_variance_all_g, num_genes, num_cells);
	cout << "data loaded \n";
      

	// Compute signal 2 noise
	double signal_to_noise_g[num_genes]; // Square root of the mean over cells of delta_gc^2/epsilon_gc^2; genes below s2n_cutoff are dropped.
	int num_genes_kept(0);
	if (s2n_cutoff > 0.0)
	  {
	    double mean_log_fold_change;
	    double log_fold_change_sample_var;
	    double mean_log_fold_change_error_sq;
	    
	    for(g=0;g<num_genes;g++)
	      {
		signal_to_noise_g[g] = 0;
		for(c=0;c<num_cells;c++)
		  {
		    signal_to_noise_g[g] += log_fold_change_all_gc[g][c]*log_fold_change_all_gc[g][c]/log_fold_change_error_sq_all_gc[g][c];
		  }
	        signal_to_noise_g[g] /= ((double) num_cells);
		signal_to_noise_g[g] = sqrt(signal_to_noise_g[g]);
		
		if(signal_to_noise_g[g] >= s2n_cutoff)
		  {
		    num_genes_kept += 1;
		  }
	      }
	  }
	else
	  {
	    num_genes_kept = num_genes;
	  }

	cout << "Done s2n determination. Number of genes left " << num_genes_kept << "\n";

	
	double **log_fold_change_gc = new double *[num_genes_kept]; // <delta_gc> of the genes above the cut-off.
	for(g=0;g<num_genes_kept;g++)
	  log_fold_change_gc[g] = new double[num_cells];

	double **log_fold_change_error_sq_gc = new double *[num_genes_kept]; // epsilon_gc^2 of the genes above the cut-off.
	for(g=0;g<num_genes_kept;g++)
	  log_fold_change_error_sq_gc[g] = new double[num_cells];

	double *gene_variance_g = new double [num_genes_kept]; // v_g of the genes above the cut-off.
	int *original_gene_index_g = new int [num_genes_kept]; // Row in the input files of each gene above the cut-off.

	if (s2n_cutoff > 0.0)
	  {
	    // Copy genes with signal to noise above cut-off
	    int k = 0;
	    for(g=0;g<num_genes;g++)
	      {
		if(signal_to_noise_g[g] >= s2n_cutoff)
		  {
		    for(c=0;c<num_cells;c++)
		      {
			log_fold_change_gc[k][c] = log_fold_change_all_gc[g][c];
			log_fold_change_error_sq_gc[k][c] = log_fold_change_error_sq_all_gc[g][c];
		      }
		    gene_variance_g[k] = gene_variance_all_g[g];
		    original_gene_index_g[k] = g;
		    k++;
		  }
	      }
	}
	else
	  {
	    for(g=0;g<num_genes;g++)
	      {
		for(c=0;c<num_cells;c++)
		  {
		    log_fold_change_gc[g][c] = log_fold_change_all_gc[g][c];
		    log_fold_change_error_sq_gc[g][c] = log_fold_change_error_sq_all_gc[g][c];
		  }
		gene_variance_g[g] = gene_variance_all_g[g];
		original_gene_index_g[g] = g;
	      }
	  }
	
	// delete temporary variables
	for(g=0;g<num_genes;g++)
	  {
	    delete[] log_fold_change_all_gc[g];
	    delete[] log_fold_change_error_sq_all_gc[g];
	  }
	delete[] log_fold_change_all_gc;
	delete[] log_fold_change_error_sq_all_gc;
	delete[] gene_variance_all_g;
	

	num_genes = num_genes_kept;
	//Initialize space for simple and full correlations
	cout << "Calculating correlations... ";
	double **corr_simple_gg = new double *[num_genes]; // Pearson correlation of <delta_gc> between every pair of genes.
	double **corr_full_gg = new double *[num_genes]; // Correlation estimated with the error bars of <delta_gc> taken into account.
	double **corr_error_gg = new double *[num_genes]; // Error bar on corr_full_gg.
	for(g=0; g<num_genes; g++)
	  {
	    corr_simple_gg[g] = new double[num_genes];
	    corr_full_gg[g] = new double[num_genes];
	    corr_error_gg[g] = new double [num_genes];
	  }


	//loop over all pairs to get simple correlations first
	int i,j;
	for(i=0;i<num_genes;i++)
	  {
	    for(j=(i+1);j<num_genes;j++)
	      {
		corr_simple_gg[i][j] = get_correlation_simple(log_fold_change_gc[i],log_fold_change_gc[j],num_cells);
		//cout << " genes " << i << " and " << j << " Rsimple " << Rs[i][j] << "\n";
	      }
	  }
	cout << "done with simple correlations \nDoing full correlations with errobars...";

       
	//Now remove the shrinkage from the prior, in place (SI Eqs. 62, 63):
	//log_fold_change_gc becomes delta~_gc and log_fold_change_error_sq_gc becomes eta_gc^2.
	double prior_removal_factor; // v_g/(v_g - epsilon_gc^2), which turns <delta_gc> into delta~_gc and epsilon_gc^2 into eta_gc^2 (SI Eqs. 62, 63).
	for(g=0;g<num_genes;g++)
	  {
	    for(c=0;c<num_cells;c++)
	      {
		if(gene_variance_g[g]-log_fold_change_error_sq_gc[g][c] > 0.0)
		  prior_removal_factor = gene_variance_g[g]/(gene_variance_g[g]-log_fold_change_error_sq_gc[g][c]);
		else
		  prior_removal_factor = gene_variance_g[g]/0.000001;
		
		log_fold_change_gc[g][c] = log_fold_change_gc[g][c]*prior_removal_factor;
		log_fold_change_error_sq_gc[g][c] = log_fold_change_error_sq_gc[g][c]*prior_removal_factor;
	      }
	  }

	//Now do the full correlation calculations
	for(i=0;i<num_genes;i++)
	  {
	    //cout << "gene number " << i << "\n";
	    for(j=(i+1);j<num_genes;j++)
	      {
		get_correlation_full(&(corr_full_gg[i][j]),&(corr_error_gg[i][j]),log_fold_change_gc[i],log_fold_change_error_sq_gc[i],gene_variance_g[i],log_fold_change_gc[j],log_fold_change_error_sq_gc[j],gene_variance_g[j],corr_simple_gg[i][j],num_cells);
	      }
	  }
	cout << "done\n";
	
	// Print correlations
	string out_file = sanity_folder + "gene_correlations";
	
	if(s2n_cutoff > 0.0)
	  {
	    string s2n_str = to_string(s2n_cutoff);
	    // Strip trailing zeros from the fractional part only, then the '.' itself.
	    // Stopping at the '.' is what keeps 10 -> "10" rather than "1".
	    if(s2n_str.find('.') != string::npos)
	      {
	        while(s2n_str.back() == '0')
	          s2n_str = s2n_str.substr(0, s2n_str.size()-1);
	        if(s2n_str.back() == '.')
	          s2n_str = s2n_str.substr(0, s2n_str.size()-1);
	      }
	    out_file += "_s2n_gt_" + s2n_str;
	  }
	out_file += ".txt";

	cout << "Printing results in " << out_file;
	print_results(out_file,corr_simple_gg,corr_full_gg,corr_error_gg,original_gene_index_g,num_genes);
  	cout << "\n";
	
	return 0;
}

//functions

double get_correlation_simple(double *log_fold_change_i_c,double *log_fold_change_j_c, int num_cells){
  int c;
  double mean_log_fold_change_i, mean_log_fold_change_j, sample_var_i, sample_var_j, sample_cov, corr; // Means, variances and covariance of <delta_gc> over the cells for the two genes; corr is their Pearson correlation.

  mean_log_fold_change_i = 0;
  mean_log_fold_change_j = 0;
  sample_var_i = 0;
  sample_var_j = 0;
  sample_cov = 0;
  for(c=0;c<num_cells;++c)
    {
      mean_log_fold_change_i += log_fold_change_i_c[c];
      mean_log_fold_change_j += log_fold_change_j_c[c];
      sample_var_i += log_fold_change_i_c[c]*log_fold_change_i_c[c];
      sample_var_j += log_fold_change_j_c[c]*log_fold_change_j_c[c];
      sample_cov += log_fold_change_i_c[c]*log_fold_change_j_c[c];
    }
  mean_log_fold_change_i /= ((double) num_cells);
  mean_log_fold_change_j /= ((double) num_cells);
  sample_var_i /= ((double) num_cells);
  sample_var_j /= ((double) num_cells);
  sample_cov /= ((double) num_cells);
  sample_var_i -= mean_log_fold_change_i*mean_log_fold_change_i;
  sample_var_j -= mean_log_fold_change_j*mean_log_fold_change_j;
  sample_cov -= mean_log_fold_change_i*mean_log_fold_change_j;
  corr = sample_cov/sqrt(sample_var_i*sample_var_j);

  return corr;
}


// Derivative of the log-likelihood of the correlation corr between two genes.
// The three arrays are defined in get_correlation_full.
double log_likelihood_slope(double corr, double *z_product_c,double *variance_inflation_product_c, double *cross_term_c, int num_cells){
  int c;
  double slope = 0.0;
  double inv_denominator;

  for(c=0;c<num_cells;++c)
    {
      inv_denominator = 1.0/(variance_inflation_product_c[c]-corr*corr);
      slope += corr*inv_denominator;
      slope += z_product_c[c]*inv_denominator;
      slope -= corr*(cross_term_c[c]-2.0*corr*z_product_c[c])*inv_denominator*inv_denominator;
    }

  return slope;
}


void get_correlation_full(double *corr_full, double *corr_error, double *log_fold_change_i_c, double *log_fold_change_error_sq_i_c, double gene_variance_i, double *log_fold_change_j_c, double *log_fold_change_error_sq_j_c, double gene_variance_j, double corr_simple, int num_cells){

  //These arrays hold simplified stats for the two genes
  static double *variance_inflation_product_c = new double[num_cells]; // variance_inflation_i * variance_inflation_j for each cell.
  static double *cross_term_c = new double[num_cells]; // z_i^2 * variance_inflation_j + z_j^2 * variance_inflation_i for each cell.
  static double *z_product_c = new double[num_cells]; // z_i * z_j for each cell.

  int c;
  double gene_sd_i = sqrt(gene_variance_i); // Standard deviation sqrt(v_g) of gene i.
  double gene_sd_j = sqrt(gene_variance_j);
  //Calculate the arrays
  for(c=0;c<num_cells;++c)
    {
      double variance_inflation_i = 1.0 + log_fold_change_error_sq_i_c[c]/gene_variance_i; // 1 + eta_gc^2/v_g: the factor by which the noise inflates the variance of delta~_gc.
      double variance_inflation_j = 1.0 + log_fold_change_error_sq_j_c[c]/gene_variance_j;
      double z_i = log_fold_change_i_c[c]/gene_sd_i; // delta~_gc in units of sqrt(v_g).
      double z_j = log_fold_change_j_c[c]/gene_sd_j;

      z_product_c[c] = z_i * z_j;
      variance_inflation_product_c[c] = variance_inflation_i * variance_inflation_j;
      cross_term_c[c] = z_i*z_i*variance_inflation_j + z_j*z_j*variance_inflation_i;
    }
  
  //Get the derivative of the log-likelihood at the simple guess corr_simple
  double slope_at_simple = log_likelihood_slope(corr_simple,z_product_c,variance_inflation_product_c,cross_term_c,num_cells);

  double corr_low, corr_high,slope_at_bound;
  if(slope_at_simple > 0)
    {
      corr_low = corr_simple;
      corr_high = 1.0;

      slope_at_bound = log_likelihood_slope(corr_high,z_product_c,variance_inflation_product_c,cross_term_c,num_cells);
      if(slope_at_bound > 0)
	{
	  //Optimum is at r=1. We estimate error-bar by by slope of log-lik at r=1.0 (i.e. assume exponential dist around r=1)
	  *corr_full = 1.0;
	  *corr_error = -log(1.0-0.159)/slope_at_bound;//note 0.156 is z-score at 1 sigma of a Gaussian. This is fairly arbitrary
	  return;
	}
    }
  else
    {
      corr_high = corr_simple;
      corr_low = -1.0;

      slope_at_bound = log_likelihood_slope(corr_low,z_product_c,variance_inflation_product_c,cross_term_c,num_cells);
      if(slope_at_bound < 0)
	{
	  *corr_full = -1.0;
	  *corr_error = log(1-0.159)/slope_at_bound;
	  return;
	}
    }

  double corr_tolerance = 0.0025;
  double slope_at_mid, corr_mid;
  while(corr_high-corr_low > corr_tolerance)
    {
      //derivative at the middle
      corr_mid = 0.5*(corr_low+corr_high);
      slope_at_mid = log_likelihood_slope(corr_mid,z_product_c,variance_inflation_product_c,cross_term_c,num_cells);

      if(slope_at_mid > 0)
	corr_low = corr_mid;
      else
	corr_high = corr_mid;
    }

  double corr_estimate = 0.5*(corr_low+corr_high);
  *corr_full = corr_estimate;

  //Now calculate the error bar (expression is a bit messy)
  double log_likelihood_curvature = 0;
  double numerator, inv_denominator;
  for(c=0;c<num_cells;++c)
    {
      numerator = variance_inflation_product_c[c]*variance_inflation_product_c[c] + 6.0 *z_product_c[c]*variance_inflation_product_c[c]*corr_estimate+2.0*z_product_c[c]*corr_estimate*corr_estimate*corr_estimate-corr_estimate*corr_estimate*corr_estimate*corr_estimate-cross_term_c[c]*(variance_inflation_product_c[c]+3.0*corr_estimate*corr_estimate);
      inv_denominator = 1.0/(variance_inflation_product_c[c]-corr_estimate*corr_estimate);
      log_likelihood_curvature += numerator*inv_denominator*inv_denominator*inv_denominator;
    }
  *corr_error = 1/sqrt(-log_likelihood_curvature);

  return;
}

void load_matrix(string sanity_folder, double **log_fold_change_gc, double **log_fold_change_error_sq_gc, double *gene_variance_g, int num_genes, int num_cells){

	string path;
	FILE *file;
	char line[3000000];
	char *token;
	char *retval;
	int g, c;
	double log_fold_change_error;

	// Read delta
	path = sanity_folder + "delta.txt";
	file = (FILE *) fopen(path.c_str(),"r");
	if(file == NULL)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",path.c_str());
	    exit(EXIT_FAILURE);
	  }
	for (g = 0; g < num_genes; g++)
	  {
	    retval = fgets(line,3000000,file);
	    if(retval == NULL)
	      {
		fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
		return;
	      }
	    token = strtok(line," \t");
	    for (c=0; c<num_cells; c++)
	      {	
		if(token != NULL)
		  {
		    log_fold_change_gc[g][c] = atof(token);
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
	fclose(file);

	// Read epsilon
	path = sanity_folder + "d_delta.txt";
	file = (FILE *) fopen(path.c_str(),"r");
	if(file == NULL)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",path.c_str());
	    exit(EXIT_FAILURE);
	  }
	for (g = 0; g < num_genes; g++)
	  {
	    retval = fgets(line,3000000,file);
	    if(retval == NULL)
	      {
		fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
		return;
	      }
	    token = strtok(line," \t");
	    for (c=0; c<num_cells; c++)
	      {	
		if(token != NULL){
		  log_fold_change_error = atof(token);
		  log_fold_change_error_sq_gc[g][c] = log_fold_change_error*log_fold_change_error;
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
	fclose(file);

	
	path = sanity_folder + "variance.txt";
	file = (FILE *) fopen(path.c_str(),"r");
	if(file == NULL)
	  {
	    fprintf(stderr,"Cannot open input file %s\n",path.c_str());
	    exit(EXIT_FAILURE);
	  }
	for (g = 0; g < num_genes; g++)
	  {
	    retval = fgets(line,3000000,file);
	    if(retval == NULL)
	      {
		fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
		return;
	      }
	    token = strtok(line," \t");
	    if(token != NULL)
	      {
		gene_variance_g[g] = atof(token);
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
	fclose(file);

	return;
}

void print_results(string out_file, double **corr_simple_gg, double **corr_full_gg, double **corr_error_gg, int *original_gene_index_g, int num_genes){

  	
	int i, j;
	string file_name(out_file);
	ofstream outfile(file_name);
	for (i=0; i<(num_genes-1); i++)
	  {
	    for (j=i+1; j<num_genes; j++)
	      {
		outfile << original_gene_index_g[i] << "\t" << original_gene_index_g[j] << "\t" << corr_simple_gg[i][j] << "\t" << corr_full_gg[i][j] << "\t" << corr_error_gg[i][j] << "\n";
	      }
	  }
	outfile.close();
}


void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, int &num_genes, int &num_cells){

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
            cout << "Sanity_gene_correlation version " << SANITY_VERSION << "\n";
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
    num_genes = -1;
    num_cells = -1;
    ifstream in(sanity_folder + "delta.txt");
    if (in.is_open()){
      num_genes = 0;
      num_cells = 0;
      string line;
      
      // get first line and count columns
      getline(in, line);
      ++num_genes;
      istringstream iss(line);
      do{
	string sub;
	iss >> sub;
	if (sub.length())
	  ++num_cells;
      }
      while(iss);
      
      // Now count genes in remaining lines
      while ( getline(in, line) ){
	++num_genes;
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

