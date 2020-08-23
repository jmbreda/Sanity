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
// g++ -std=c++11 -O2 -ffast-math -O3 -fopenmp compute_distance.cpp -o Sanity_distance

// Distances functions :
// with error bars:
void get_distance_errorbar(double **delta, double **epsilon2, double *variance, double **D, int N_gene, int N_cell, int N_threads);
void get_Di_errorbar(double **delta, double **epsilon2, double **alpha_v, double *D, int N_gene, int N_cell, int i, int N_bin);
// without errorbars:
void get_distance_euclidean(double **delta, double **D, int N_gene, int N_cell, int N_threads);
void get_Di_euclidean(double **delta, double *D, int N_gene, int N_cell, int i);

// I/O
void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, bool &with_error_bar, int &N_threads, int &N_gene, int &N_cell);
static void show_usage(void);
void load_matrix(string sanity_folder, double **delta, double **epsilon2, double *variance, int N_gene, int N_cell);
void print_distance(string out_file, double **D, int N_cell);

int main(int argc, char** argv){


	// Declare arguments :
	string sanity_folder("");
	double s2n_cutoff = 1.0;
	bool with_error_bar(true);
	int N_threads(4);
	int N_gene(0);
	int N_cell(0);
	parse_argv(argc, argv, sanity_folder, s2n_cutoff, with_error_bar, N_threads, N_gene, N_cell);

	cout << "Sanity folder :\n\t" << sanity_folder << "\n";
	cout << "Signal/noise cut-off :" << s2n_cutoff << "\n";
	if (with_error_bar)
		cout << "With error bars";
	else
		cout << "Without error bars";
	cout << "\n";
	cout << "Nr. of gene :" << N_gene << "\n";
	cout << "Nr. of cell :" << N_cell << "\n";
	cout << "Using " << N_threads << " cpus\n";


	// Declare temporary variables
	int g, c;
	double **delta_all = new double *[N_cell];
	for(c=0;c<N_cell;c++)
		delta_all[c] = new double[N_gene];

	double **epsilon2_all = new double *[N_cell];
	for(c=0;c<N_cell;c++)
		epsilon2_all[c] = new double[N_gene];
	
	double *variance_all = new double [N_gene];

	// Get data
	cout << "Get data... ";
	load_matrix(sanity_folder, delta_all, epsilon2_all, variance_all, N_gene, N_cell);	

	// Compute signal 2 noise
	bool idx_gene[N_gene];
	int N_gene_cutoff(0);
	if (s2n_cutoff > 0.0){
	
		double mean_delta;
		double var_delta;
		double mean_epsilon2;
		for(g=0;g<N_gene;g++){
			
			mean_delta = 0.0;
			for(c=0;c<N_cell;c++){
				mean_delta += delta_all[c][g];
			}
			mean_delta /= (double)N_cell;

			var_delta = 0.0;
			for(c=0;c<N_cell;c++){
				var_delta += (delta_all[c][g] - mean_delta)*(delta_all[c][g] - mean_delta);
			}
			var_delta /= (double)N_cell - 1.0 ;
	
			mean_epsilon2 = 0.0;
			for(c=0;c<N_cell;c++){
				mean_epsilon2 += epsilon2_all[c][g];
			}
			mean_epsilon2 /= (double)N_cell;

			if(var_delta/mean_epsilon2 >= s2n_cutoff){
				idx_gene[g] = true;
				N_gene_cutoff += 1;
			}
			else{
				idx_gene[g] = false;
			}
		}
	}else{
		N_gene_cutoff = N_gene;
	}
	
	double **delta = new double *[N_cell];
	for(c=0;c<N_cell;c++)
		delta[c] = new double[N_gene_cutoff];

	double **epsilon2 = new double *[N_cell];
	for(c=0;c<N_cell;c++)
		epsilon2[c] = new double[N_gene_cutoff];

	double *variance = new double [N_gene_cutoff];

	if (s2n_cutoff > 0.0){
		// Copy genes with signal to noise above cut-off
		int k = 0;
		for(g=0;g<N_gene;g++){
			if(idx_gene[g]){
				for(c=0;c<N_cell;c++){
					delta[c][k] = delta_all[c][g];
					epsilon2[c][k] = epsilon2_all[c][g];
				}
				variance[k] = variance_all[g];
				k++;
			}
		}
	}else{
		for(g=0;g<N_gene;g++){
            for(c=0;c<N_cell;c++){
                delta[c][g] = delta_all[c][g];
                epsilon2[c][g] = epsilon2_all[c][g];
            }
            variance[g] = variance_all[g];
        }
	}

	N_gene = N_gene_cutoff;
	cout << "N_gene after cut-off :" << N_gene << "\n"; 

	// delete temporary variables
	for(c=0;c<N_cell;c++){
		delete[] delta_all[c];
		delete[] epsilon2_all[c];
	}
	delete[] delta_all;
	delete[] epsilon2_all;
	delete[] variance_all;

	// Get distance
	cout << "Get distance... ";
	double **D = new double *[N_cell];
	for(c=0; c<N_cell; c++)
		D[c] = new double[N_cell];

	int i,j;
	for(i=0;i<N_cell;i++){
		for(j=0;j<N_cell;j++){
			D[i][j] = 0.0;
		}
	}
	if(with_error_bar)
		get_distance_errorbar(delta,epsilon2,variance,D,N_gene,N_cell,N_threads);
	else
		get_distance_euclidean(delta,D,N_gene,N_cell,N_threads);
  	cout << "\n";

	// Print distance
	string out_file = sanity_folder + "cell_cell_distance";
	if(with_error_bar)
		out_file += "_with_errorbar";
	else
		out_file += "_euclidean";
	if(s2n_cutoff > 0.0){
		string s2n_str = to_string(s2n_cutoff);
		while(s2n_str.back() == '0' || s2n_str.back()=='.')
			s2n_str = s2n_str.substr(0, s2n_str.size()-1);
		out_file += "_s2n_gt_" + s2n_str;
	}
	out_file += ".txt";

	cout << "Print distance in " << out_file;
	print_distance(out_file, D, N_cell);
  	cout << "\n";

	return 0;
}

void get_distance_errorbar(double **delta, double **epsilon2, double *variance, double **D, int N_gene, int N_cell, int N_threads){

	// Rescaling of delta and epsilon2
	int g, c;
	double factor;
	for(g=0;g<N_gene;g++){
		for(c=0;c<N_cell;c++){
			if( variance[g]-epsilon2[c][g] > 0.0 )
				factor = variance[g]/(variance[g]-epsilon2[c][g]);
			else
				factor = variance[g]/0.000001;

			delta[c][g] = delta[c][g]*factor;
			epsilon2[c][g] = epsilon2[c][g]*factor;
		}
	}

	// define variance_gene*alpha_bin*Alpha_gene
	int k;
	int N_bin(401);
	double da(0.005);
	double a;

	double **alpha_v = new double *[N_gene];
	for(g=0;g<N_gene;g++)
		alpha_v[g] = new double [N_bin];

	for(g=0;g<N_gene;g++){
		a = 0.0;
		for(k=0;k<N_bin;k++){
			alpha_v[g][k] = a*variance[g];
			a += da;
		}
	}

	int i, N_parloop;
	if ( (N_cell % 2) == 0 ){
		N_parloop = N_cell/2;
	}else{
		N_parloop =(N_cell-1)/2;
		i = N_parloop;
		get_Di_errorbar(delta,epsilon2,alpha_v,D[i],N_gene,N_cell,i,N_bin);
	}
	#pragma omp parallel for num_threads(N_threads)
	for (i=0; i<N_parloop; i++){
		get_Di_errorbar(delta,epsilon2,alpha_v,D[i],N_gene,N_cell,i,N_bin);
		get_Di_errorbar(delta,epsilon2,alpha_v,D[N_cell-1-i],N_gene,N_cell,N_cell-i-1,N_bin);
	}
	return;
}

void get_Di_errorbar(double **delta, double **epsilon2, double **alpha_v, double *D, int N_gene, int N_cell, int i, int N_bin){
	
	int j, g, k;
	double q, f;
	double x2[N_gene];
	double eps2[N_gene];
	double lik[N_bin];
	double lik_max, lik_tot;
	for (j=i+1; j<N_cell; j++){
		
		// Compute 
		// x_g^2 = (delta_gi - delta_gj)^2
		// eps2 = epsilon2_g1 + epsilon2_g2
		// f = alpha_v/(eps2 + alpha_v)
		for(g=0;g<N_gene;g++){
			x2[g] = (delta[i][g] - delta[j][g])*(delta[i][g] - delta[j][g]);
			eps2[g] = epsilon2[i][g] + epsilon2[j][g];
		}

		// get likelihood
		lik_max = -1000000000.0;
		for(k=0;k<N_bin;k++){
			// log-likelihood over all genes
			lik[k] = 0.0;
			for(g=0;g<N_gene;g++){
				q= eps2[g] + alpha_v[g][k];
				lik[k] -= .5*x2[g]/q;
				lik[k] -= .5*log(q);
			}
			// keep  track of the maximum log-likelihood
			if(lik[k] > lik_max)
				lik_max = lik[k];
		}

		// normalize
		lik_tot = 0.0;
		for(k=0;k<N_bin;k++){
			lik[k] -= lik_max;
			// turn into likelihood
			lik[k] = exp(lik[k]);
			// keep track of normalization
			lik_tot += lik[k];
		}

		// calculate expectation values
		D[j] = 0.0;
		for(k=0;k<N_bin;k++){
			// posterior of this alpha
			lik[k] /= lik_tot;
			for(g=0;g<N_gene;g++){
				// This is f_g(alpha)
				f = alpha_v[g][k]/(alpha_v[g][k] + eps2[g]);
				// this is contribution to distance squared
				D[j] += lik[k]*(f*f*x2[g] + f*eps2[g]);
			}
		}
		D[j] = sqrt(D[j]);
	}
	return;
}


void get_distance_euclidean(double **delta, double **D, int N_gene, int N_cell, int N_threads){

    int i, N_parloop;
    if ( (N_cell % 2) == 0 ){
        cout << "N_cell even\n";
        N_parloop = N_cell/2;
    }else{
        cout << "N_cell odd\n";
        N_parloop =(N_cell-1)/2;
        i = N_parloop;
        get_Di_euclidean(delta,D[i],N_gene,N_cell,i);
    }
    #pragma omp parallel for num_threads(N_threads)
    for (i=0; i<N_parloop; i++){
        get_Di_euclidean(delta,D[i],N_gene,N_cell,i);
        get_Di_euclidean(delta,D[N_cell-1-i],N_gene,N_cell,N_cell-i-1);
    }
    return;
}

void get_Di_euclidean(double **delta, double *D, int N_gene, int N_cell, int i){

    int j, g;
    double diff;
    for (j=i+1; j<N_cell; j++){
        D[j] = 0.0;
        for (g=0; g<N_gene; g++){
            diff = delta[i][g] - delta[j][g];
            D[j] += diff*diff;
        }
        D[j] = sqrt(D[j]);
    }

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
	if(infp == NULL){
		fprintf(stderr,"Cannot open input file %s\n",my_file.c_str());
		exit(EXIT_FAILURE);
	}
	for (g = 0; g < N_gene; g++){
		retval = fgets(ss,3000000,infp);
		if(retval == NULL){
			fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
			return;
		}
		token = strtok(ss," \t");
		for (c=0; c<N_cell; c++){	
			if(token != NULL){
				delta[c][g] = atof(token);
			}else{
				fprintf(stderr,"Error: not enough fields on line number %d:\n",g);
				return; 
			}
			token = strtok(NULL," \t");
		}
		token = strtok(NULL," \t");
		if(token != NULL){
			fprintf(stderr,"Error: too many fields on line number %d:\n",g);
		}
	}
	fclose(infp);

	// Read epsilon
	my_file = sanity_folder + "d_delta.txt";
	infp = (FILE *) fopen(my_file.c_str(),"r");
	if(infp == NULL){
		fprintf(stderr,"Cannot open input file %s\n",my_file.c_str());
		exit(EXIT_FAILURE);
	}
	for (g = 0; g < N_gene; g++){
		retval = fgets(ss,3000000,infp);
		if(retval == NULL){
			fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
			return;
		}
		token = strtok(ss," \t");
		for (c=0; c<N_cell; c++){	
			if(token != NULL){
				epsilon = atof(token);
				epsilon2[c][g] = epsilon*epsilon;
			}else{
				fprintf(stderr,"Error: not enough fields on line number %d:\n",g);
				return; 
			}
			token = strtok(NULL," \t");
		}
		token = strtok(NULL," \t");
		if(token != NULL){
			fprintf(stderr,"Error: too many fields on line number %d:\n",g);
		}
	}
	fclose(infp);

	// Read variance
	my_file = sanity_folder + "variance.txt";
	infp = (FILE *) fopen(my_file.c_str(),"r");
	if(infp == NULL){
		fprintf(stderr,"Cannot open input file %s\n",my_file.c_str());
		exit(EXIT_FAILURE);
	}
	for (g = 0; g < N_gene; g++){
		retval = fgets(ss,3000000,infp);
		if(retval == NULL){
			fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
			return;
		}
		token = strtok(ss," \t");
		if(token != NULL){
			variance[g] = atof(token);
		}else{
			fprintf(stderr,"Error: not enough fields on line number %d:\n",g);
			return;
		}
		token = strtok(NULL," \t");
		if(token != NULL){
			fprintf(stderr,"Error: too many fields on line number %d:\n",g);
		}
	}
	fclose(infp);

	return;
}


void print_distance(string out_file, double **D, int N_cell){
	
	int i, j;
	string file_name(out_file);
	ofstream outfile(file_name);
	for (i=0; i<(N_cell-1); i++){
		for (j=i+1; j<N_cell; j++)
			outfile << D[i][j] << "\n";
	}
	outfile.close();
}


void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, bool &with_error_bar, int &N_threads, int &N_gene, int &N_cell){

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
    for(i=1;i<argc;i++){
        if (argv[i] == get_version[0] || argv[i] == get_version[1]){
            cout << "v1.0" << "\n";
            exit(0);
        }
    }

    int N_param(4);
    string error_bar("true");
    string to_find[4][2] = {{"-f", "--folder"},
                            {"-s2n", "--signal_to_noise_cutoff"},
                            {"-err","--with_error_bars"},
                            {"-n", "--n_threads"}};

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
				if(j==2) error_bar = argv[idx+1];
                if(j==3) N_threads = atoi(argv[idx+1]);

                // add '/' to out_folder if not already
                if( j == 0 && sanity_folder.back() != '/' )
                    sanity_folder = sanity_folder + '/';
            }
        }
        if (idx == 0 && j == 0){
            cerr << "Error in argument parsing :\n"
            << "missing input file name\n";
            show_usage();
        }
    }
    
	if ( error_bar == "false" || error_bar == "0" )
        with_error_bar = false;

    // Get number of gene and number of cells
	string command;
	int out;
	string file_name = sanity_folder + "n_gene_n_cell";
    command = "wc -l " + sanity_folder + "delta.txt|cut -f1 -d' ' > " + file_name;
    out = system(command.c_str());
	command = "head -1 " + sanity_folder + "delta.txt|wc -w >> " + file_name;
    out = system(command.c_str());
	ifstream myfile (file_name);
	string line;
    if (myfile.is_open()){
        getline (myfile,line);
        N_gene = atoi(line.c_str());
        getline (myfile,line);
        N_cell = atoi(line.c_str());
        myfile.close();
        command = "rm " + file_name;
        out = system(command.c_str());
    }
    else{
        cerr << "Unable to open " + file_name + "\n";
    }
}

static void show_usage(void)
{
    cerr << "Usage: Sanity_distance <option(s)> SOURCES\n"
         << "Options:\n"
         << "\t-h,--help\t\tShow this help message\n"
         << "\t-v,--version\t\tShow the current version\n"
         << "\t-f,--folder\t\tSpecify the input folder with extended output from Sanity\n"
         << "\t-s2n,--signal_to_noise_cutoff\tMinimal signal/noise of genes to include in the distance calculation (default: 1.0)\n"
         << "\t-err,--with_error_bars\tCompute cell-cell distance taking the errobar epsilon into account (default: 1)\n"
         << "\t-n,--n_threads\t\tSpecify the number of threads to be used (default: 4)\n";
    exit(0);
}



