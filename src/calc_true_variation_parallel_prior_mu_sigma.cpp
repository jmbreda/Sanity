#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <fstream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <sstream>

#include <sys/stat.h>  // For mkdir on Unix-like systems
#ifdef _WIN32
    #include <direct.h>  // For _mkdir on Windows
#endif

#include <ReadInputFiles.h>
#include <FitFrac.h>
#include <Digamma_Trigamma.h>

string VERSION("1.1");

using namespace std;

struct RowComputation {
	double mu;
	double var_mu;
	std::vector<double> delta;
	std::vector<double> var_delta;
	double var_gene;
	std::vector<double> lik;
	double var_gene_v_ml;
	double mu_v_ml;
	double var_mu_v_ml;
	std::vector<double> delta_v_ml;
	std::vector<double> var_delta_v_ml;
};

/***Function declarations ****/
RowComputation get_gene_expression_level(const vector<double>& n_c, const vector<double>& N_c, double n, double vmin, double vmax, int C, int numbin, double a, double b, bool max_v_output, bool post_v_output);
double get_epsilon_2(double &d, double &v, double &n, double &f, double &a);
void parse_argv(int argc,char** argv, string &in_file, string &gene_name_file, string &cell_name_file, string &in_file_extension, string &out_folder, int &N_threads, bool &print_extended_output, double &vmin, double &vmax, int &numbin, int &N_charm, bool &no_norm, bool &max_v_output, bool &post_v_output);
static void show_usage(void);
std::vector<double> fetch_row(int g, const std::string& in_file, const std::string& in_file_extension, const std::vector<RowBlock>& mtx_rows, const std::vector<std::streampos>& tsv_offsets, const int& C);

// Function to print a timestamped message to stderr
void logging_debug(const std::string& msg) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::tm buf;
#ifdef _WIN32
    localtime_s(&buf, &in_time_t);
#else
    localtime_r(&in_time_t, &buf);
#endif
    std::cerr << "[" << std::put_time(&buf, "%Y-%m-%d %H:%M:%S");
    std::cerr << ',' << std::setw(3) << std::setfill('0') << ms.count() << "] ";
    std::cerr << msg << std::endl;
}

int main (int argc, char** argv){
	string in_file("");
	string gene_name_file("none");
	string cell_name_file("none");
	string in_file_extension("");
	string out_folder("");
	int N_threads(4);
	bool print_extended_output(false);
	double vmin = 0.001;
    double vmax = 50.0;
	int numbin = 160;
	int N_char;
	bool no_norm(false);
	bool max_v_output(false);
	bool post_v_output(true);
	parse_argv(argc, argv, in_file, gene_name_file, cell_name_file, in_file_extension, out_folder, N_threads, print_extended_output, vmin, vmax, numbin, N_char, no_norm, max_v_output, post_v_output);

	// count Number of genes and cells
	int G, C;
	// Number of rows in file
	int N_rows;
	// Gene idx map for mtx
	map<int,int> gene_idx;
	vector<RowBlock> mtx_rows;
	vector<std::streampos> tsv_offsets;
	vector<double> N_c; // total counts per cell
	vector<double> n; // total counts per gene
	vector<string> gene_names;
	vector<string> cell_names;

	if (in_file_extension == "mtx"){
		Get_G_C_MTX(in_file, N_rows, G, C, gene_idx, mtx_rows, N_c, n);
		gene_names = Read_GeneNames(gene_name_file, gene_idx, G);
		cell_names = Read_CellNames(cell_name_file);
		logging_debug("There were " + std::to_string(N_rows) + " rows");
	} else {
		Get_G_C_UMIcountMatrix(in_file, N_rows, G, C, N_char, tsv_offsets, N_c, n, cell_names, gene_names);
		logging_debug("There were " + std::to_string(N_rows) + " rows");
	}
	logging_debug("There were " + std::to_string(G) + " genes and " + std::to_string(C) + " cells");

	int g, c, k;

	// Remove the total UMI correction if no cell size normalization option is true
	if (no_norm){
		logging_debug("No cell size normalization will be performed");

		// get mean count per cell
		double mean_N_c = 0;
		for(c=0;c<C;++c){
			mean_N_c += N_c[c];
		}
		mean_N_c /= C;

		// Now replce N_c by N
		for(c=0;c<C;++c){
			N_c[c] = mean_N_c;
		}
	}

	// alpha and beta of gamma prior on mu
	double a = 1.0;
	double b = 0.0;
	double deltav = log(vmax/vmin)/((double) numbin-1);

	/*
	for (int g = 0; g < G; ++g) {
			std::vector<double> n_c_g = fetch_row(g, in_file, in_file_extension, mtx_rows, tsv_offsets, C);
			RowComputation result = get_gene_expression_level(n_c_g, N_c, n[g], vmin, vmax, C, numbin, a, b, max_v_output, post_v_output);
			cout << g << '\n';
	}
	*/
	
	//open files for writing
	std::ofstream out_exp_lev, out_d_exp_lev, out_mu, out_dmu, out_var_gene, \
		out_delta, out_ddelta, out_lik, out_exp_lev_v_ml, out_d_exp_lev_v_ml, \
		out_mu_v_ml, out_dmu_v_ml, out_var_gene_v_ml, out_delta_v_ml, out_ddelta_v_ml,
		out_gene, out_cell;

	if (post_v_output) {
		std::ofstream(out_folder + "log_transcription_quotients.txt", std::ios::trunc).close();
		std::ofstream(out_folder + "ltq_error_bars.txt", std::ios::trunc).close();

		out_exp_lev.open(out_folder + "log_transcription_quotients.txt", std::ios::app);
		out_d_exp_lev.open(out_folder + "ltq_error_bars.txt", std::ios::app);
		
		out_exp_lev << "GeneID";
		out_d_exp_lev << "GeneID";
		for(c=0;c<C;c++){
            out_exp_lev << "\t" << cell_names[c].c_str();
            out_d_exp_lev << "\t" << cell_names[c].c_str();
        }
        out_exp_lev << "\n";
        out_d_exp_lev << "\n";
		if ( print_extended_output ){
			std::ofstream(out_folder + "geneID.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "cellID.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "mu.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "d_mu.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "variance.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "delta.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "d_delta.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "likelihood.txt", std::ios::trunc).close();

			out_gene.open(out_folder + "geneID.txt", std::ios::app);
			out_cell.open(out_folder + "cellID.txt", std::ios::app);
			out_mu.open(out_folder + "mu.txt", std::ios::app);
			out_dmu.open(out_folder + "d_mu.txt", std::ios::app);
			out_var_gene.open(out_folder + "variance.txt", std::ios::app);
			out_delta.open(out_folder + "delta.txt", std::ios::app);
			out_ddelta.open(out_folder + "d_delta.txt", std::ios::app);
			out_lik.open(out_folder + "likelihood.txt", std::ios::app);

			out_lik << "Variance";
			for(k=0;k<(numbin);++k){
                out_lik << "\t" << vmin*exp(deltav*k);
            }
            out_lik << "\n";
		}
	}
	if (max_v_output) {
		std::ofstream(out_folder + "log_transcription_quotients_vmax.txt", std::ios::trunc).close();
		std::ofstream(out_folder + "ltq_error_bars_vmax.txt", std::ios::trunc).close();

		out_exp_lev_v_ml.open(out_folder + "log_transcription_quotients_vmax.txt", std::ios::app);
		out_d_exp_lev_v_ml.open(out_folder + "ltq_error_bars_vmax.txt", std::ios::app);
		
		out_exp_lev_v_ml << "GeneID";
		out_d_exp_lev_v_ml << "GeneID";
		for(c=0;c<C;c++){
            out_exp_lev_v_ml << "\t" << cell_names[c].c_str();
            out_d_exp_lev_v_ml << "\t" << cell_names[c].c_str();
        }
        out_exp_lev_v_ml << "\n";
        out_d_exp_lev_v_ml << "\n";
		if ( print_extended_output ){
			std::ofstream(out_folder + "geneID.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "cellID.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "mu_vmax.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "d_mu_vmax.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "variance_vmax.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "delta_vmax.txt", std::ios::trunc).close();
			std::ofstream(out_folder + "d_delta_vmax.txt", std::ios::trunc).close();

			out_gene.open(out_folder + "geneID.txt", std::ios::app);
			out_cell.open(out_folder + "cellID.txt", std::ios::app);
			out_mu_v_ml.open(out_folder + "mu_vmax.txt", std::ios::app);
			out_dmu_v_ml.open(out_folder + "d_mu_vmax.txt", std::ios::app);
			out_var_gene_v_ml.open(out_folder + "variance_vmax.txt", std::ios::app);
			out_delta_v_ml.open(out_folder + "delta_vmax.txt", std::ios::app);
			out_ddelta_v_ml.open(out_folder + "d_delta_vmax.txt", std::ios::app);
		}	
	}

	logging_debug("Fit gene expression levels");
	const clock_t begin = clock();
	#pragma omp parallel num_threads(N_threads)
	{
		#pragma omp for schedule(dynamic) ordered
		for (int g = 0; g < G; ++g) {
			std::vector<double> n_c_g = fetch_row(g, in_file, in_file_extension, mtx_rows, tsv_offsets, C);
			RowComputation result = get_gene_expression_level(n_c_g, N_c, n[g], vmin, vmax, C, numbin, a, b, max_v_output, post_v_output);
			#pragma omp ordered
			{
                          if (g == (3 * N_threads - 1)) {
			                        clock_t end = clock();
                  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                  	double estimated_running_time = elapsed_secs/(3 * N_threads) *(G-3 * N_threads)/N_threads;
 
                  	// output esimated running time
                  	string estimated_running_time_str = "Estimated running time: ";
                  	if ( estimated_running_time > 86400 ){
						std::stringstream stream;
						stream << std::fixed << std::setprecision(2) << estimated_running_time/(3600*24);
                        estimated_running_time_str += stream.str() + " days ";
                  	}else if ( estimated_running_time > 3600 ){
						std::stringstream stream;
						stream << std::fixed << std::setprecision(2) << estimated_running_time/3600;
                        estimated_running_time_str += stream.str() + " hours ";
                  	}else if ( estimated_running_time > 60 ){
						std::stringstream stream;
						stream << std::fixed << std::setprecision(2) << estimated_running_time/60;
                        estimated_running_time_str += stream.str() + " minutes ";
                  	}else{
						std::stringstream stream;
						stream << std::fixed << std::setprecision(2) << estimated_running_time;
                        estimated_running_time_str += stream.str() + " seconds ";
                  	}
					logging_debug(estimated_running_time_str);
				}
				if (g % 1000 == 0 && g > 0){
                    logging_debug("Finished " \
						+ std::to_string(g) + " genes out of " + std::to_string(G));
                  }
				if (post_v_output) {
					out_exp_lev << gene_names[g];
					out_d_exp_lev << gene_names[g];
					for (c=0;c<C;c++){
						out_exp_lev << "\t" << result.mu + result.delta[c];
						out_d_exp_lev << "\t" << sqrt(result.var_mu + result.var_delta[c]);
						if ( print_extended_output ){
							out_delta << std::fixed << std::setprecision(6) << result.delta[c];
							out_ddelta << std::fixed << std::setprecision(6) << sqrt(result.var_delta[c]);
							if ( c < C-1 ) {
								out_delta << "\t";
								out_ddelta << "\t";	
							}
						}
					}
					out_exp_lev << "\n";
					out_d_exp_lev << "\n";
					if ( print_extended_output ){
						out_delta << "\n";
						out_ddelta << "\n";
						// Write gene names
						out_gene << gene_names[g].c_str() << "\n";
						//print best fit to file : mu, delta
						// Print diagonal of invM : variance of mu, delta
						out_mu << std::fixed << std::setprecision(6) << result.mu << "\n";
						out_dmu << std::fixed << std::setprecision(6) << sqrt(result.var_mu) << "\n";
						out_var_gene << std::fixed << std::setprecision(6) << result.var_gene << "\n";
						// Write likelihood
						out_lik << gene_names[g];
						for(k=0;k<numbin;++k){
							out_lik << "\t" << result.lik[k];
						}
						out_lik << "\n";
					}
				}
				if (max_v_output) {
					out_exp_lev_v_ml << gene_names[g];
					out_d_exp_lev_v_ml << gene_names[g];
					for (c=0;c<C;c++){
						out_exp_lev_v_ml << "\t" << result.mu_v_ml + result.delta_v_ml[c];
						out_d_exp_lev_v_ml << "\t" << sqrt(result.var_mu_v_ml + result.var_delta_v_ml[c]);
						if ( print_extended_output ){
							out_delta_v_ml << std::fixed << std::setprecision(6) << result.delta_v_ml[c];
							out_ddelta_v_ml << std::fixed << std::setprecision(6) << sqrt(result.var_delta_v_ml[c]);
							if ( c < C-1 ) {
								out_delta_v_ml << "\t";
								out_ddelta_v_ml << "\t";	
							}	
						}
					}
					out_exp_lev_v_ml << "\n";
					out_d_exp_lev_v_ml << "\n";
					if ( print_extended_output ){
						out_delta_v_ml << "\n";
						out_ddelta_v_ml << "\n";
						// Write gene names
						out_gene << gene_names[g].c_str() << "\n";
						//print best fit to file : mu, delta
						// Print diagonal of invM : variance of mu, delta
						out_mu_v_ml << std::fixed << std::setprecision(6) << result.mu_v_ml << "\n";
						out_dmu_v_ml << std::fixed << std::setprecision(6) << sqrt(result.var_mu_v_ml) << "\n";
						out_var_gene_v_ml << std::fixed << std::setprecision(6) << result.var_gene_v_ml << "\n";
					}
				}
			}
		}
	}

	// save cell names
	for(c=0;c<C;c++){
		out_cell << cell_names[c].c_str() << "\n";
	}

	logging_debug("Finished fitting all genes");

    return 0;
}

RowComputation get_gene_expression_level(const vector<double>& n_c, const vector<double>& N_c, \
	double n, double vmin, double vmax, int C, int numbin, double a, double b, \
	 bool max_v_output, bool post_v_output) {
	int i, k;
    double beta,L,ldet,q,delsq,inv_v;
	double *f = new double[C];
	double **delta_v = new double *[numbin];
    double **sig2_delta_v = new double *[numbin];
	std::vector<double> lik(numbin, -1.0);
	std::vector<double> delta(C, 0.0), var_delta(C, 0.0), delta_v_ml(C, 0.0), var_delta_v_ml(C, 0.0);
	double var_mu;
	double var_gene_v_ml;
	double mu_v_ml;
	double var_mu_v_ml;

	for(k=0;k<numbin;++k){
        delta_v[k] = new double [C];
        sig2_delta_v[k] = new double [C];
    }

	/*** To compute var of delta ***/
    double *sig2_delta_c = new double [C];
    double *sig2_delta_num = new double [C];
    double *sig2_delta_den2 = new double [C];
    double sig2_delta_den1;

	double *mu_v = new double[numbin];
	double Lmax = -1e+100;
	int Lmax_ind = 0;
    double v;
	double deltav;
    deltav = log(vmax/vmin)/((double) numbin-1);

	for(k=0;k<numbin;++k){
		v = vmin * exp(deltav*k);
		beta = 1.0/((n+a)*v);
		//#pragma omp critical
		q = fitfrac(f,n_c,n,v,C,N_c,a,b);

		mu_v[k] = Psi_0(n+a)-q; /*** equation (85) ***/

		double tot = 0;
		for(i=0;i<C;++i){
			tot += f[i];
		}
		delsq = 0;
		L = -0.5*((double) C)*log(v); // (56) 1st term
		for(i=0;i<C;++i){
			delta_v[k][i] = log(f[i]) - log(N_c[i]) + q; // equation (67)
			L += n_c[i]*delta_v[k][i];// Bug fix: remove a term as in Equation 19 of Sanity paper SI
			delsq += delta_v[k][i]*delta_v[k][i];
		}
		L -= delsq/(2*v);// (56) 2nd term
		L -= (n+a)*q;//4th term in equation (56)

		// get the determinant of the matrix
		ldet = 0.0;
		for(i=0;i<C;++i){
			ldet += (f[i]*f[i])/(f[i]+beta);
		}

		ldet = log(1 - ldet);
		for(i=0;i<C;++i){
			ldet += log(f[i]+beta);
		}
		L -= 0.5*ldet;
		// substract prior with a = 1, b = 1 ( log(v^a*exp(-b*v)) = alog(v) - bv
		lik[k] = L;
		
		if(L > Lmax){
			Lmax = L;
			Lmax_ind = k;
		}

		/* compute nf^2/(nf+1/sigma^2) for each c */
		inv_v = 1.0/v;
		for(i=0;i<C;++i){
			sig2_delta_c[i] = (n+a)*f[i]*f[i]/((n+a)*f[i] + inv_v);
		}
		/* Compute the full sum of the denominator in Delta_delta and the second tern in the denominator*/
		sig2_delta_den1 = 1.0;
		for(i=0;i<C;++i){
			sig2_delta_den1 -= sig2_delta_c[i];
			sig2_delta_den2[i] = (n+a)*f[i] + inv_v;
		}
		/* compute the different terms in the numerator : remove the \tilde{c} terms */
		for(i=0;i<C;i++){
			sig2_delta_num[i] = sig2_delta_den1 + sig2_delta_c[i];
		}
		/* compute sig2_delta */
		for(i=0;i<C;++i){
			sig2_delta_v[k][i] = sig2_delta_num[i]/(sig2_delta_den1*sig2_delta_den2[i]);
		}

        // fix computation of asymmetric sig2_delta for zero count
        for(i=0;i<C;++i){
            if(n_c[i]<=0.5){
                sig2_delta_v[k][i] = get_epsilon_2(delta_v[k][i],v,n,f[i],a);
            }
        }
	}// end v bins loop

	// get normalized likelihood from loglikelihood
	double sum_L;
	sum_L = 0.0;
	for(k=0;k<numbin;k++){
		lik[k] -= Lmax;
		lik[k] = exp(lik[k]);
		sum_L += lik[k];
	}
	for(k=0;k<numbin;k++){
		lik[k] /= sum_L;
	}

	/*
	// Multiply by prior
	for(k=0;k<numbin;k++){
		v = vmin * exp(deltav*k);
		lik[k] *= v*exp(-v);
	}
	// Renormalize
	sum_L = 0;
	for(k=0;k<numbin;k++){
		sum_L += lik[k];
	}
	for(k=0;k<numbin;k++){
		lik[k] /= sum_L;
	}
	*/

	// Average delta, and mu
	double mu = 0.0;
	for(k=0;k<numbin;k++){
		mu += lik[k]*mu_v[k];
	}

	// Compute var_delta = < (mu - <mu>)^2 > + <d_mu>
	var_mu = Psi_1((double) n + 1);
	for(k=0;k<numbin;k++){
		var_mu += lik[k]*(mu_v[k] - mu)*(mu_v[k] - mu);
	}

	if (post_v_output){
		// Compute <delta> = int p(v)*delta(v) dv
		for(i=0;i<C;i++){
			delta[i] = 0.0;
			for(k=0;k<numbin;k++){
				delta[i] += lik[k]*delta_v[k][i];
			}
		}
		// Compute var_delta = < (delta - <delta>)^2 > + <sig2_delta^2>
		for(i=0;i<C;i++){
			var_delta[i] = 0.0;
			for(k=0;k<numbin;k++){
				var_delta[i] += lik[k]*(delta_v[k][i]-delta[i])*(delta_v[k][i]-delta[i]) + lik[k]*sig2_delta_v[k][i];
			}
		}
	}

	if (max_v_output){
		// Store the gene-variance that maximizes the likelihood:
		var_gene_v_ml = vmin * exp(deltav*Lmax_ind);
		// And then also the corresponding values for the LTQs etc.
		mu_v_ml = mu_v[Lmax_ind];
		var_mu_v_ml = Psi_1((double) n + 1);
		for(i=0;i<C;i++){
			delta_v_ml[i] = delta_v[Lmax_ind][i];
		}
		for(i=0;i<C;i++){
			var_delta_v_ml[i] = sig2_delta_v[Lmax_ind][i];
		}
	}

	delete[] f;
	for(k=0;k<numbin;++k){
        delete[] delta_v[k];
        delete[] sig2_delta_v[k];
    }
	delete[] delta_v;
	delete[] sig2_delta_v;
	delete[] sig2_delta_c;
    delete[] sig2_delta_num;
    delete[] sig2_delta_den2;
	delete[] mu_v;

	RowComputation result;
	result.mu = mu;
	result.var_mu = var_mu;
	if (post_v_output){
		result.delta = delta;
		result.var_delta = var_delta;
		result.var_gene = 0.0;
 		for(k=0;k<numbin;++k){
			v = vmin * exp(deltav*k);
			result.var_gene += v*lik[k];
		}
		result.lik = lik;
	}
	if (max_v_output){
		result.var_gene_v_ml = var_gene_v_ml;
		result.mu_v_ml = mu_v_ml;
		result.var_mu_v_ml = var_mu_v_ml;
		result.delta_v_ml = delta_v_ml;
		result.var_delta_v_ml = var_delta_v_ml;
	}
	return result;
}

double get_epsilon_2(double &d, double &v, double &n, double &f, double &a){

    double e;
    double dL;
    double e_low = 0.0;
    double e_high = 0.0;
    double vnf = v*(n+a)*f;
    e_high = ( -(d+vnf) + sqrt( (d+vnf)*(d+vnf) + v*(1.0 + vnf) ) )/(1.0 + vnf);

    // bisection method :
    double tol = 0.0000001;
    double diff = 1.0;
    while(diff > tol){
        e = (e_high+e_low)/2.0;
        dL = e*(2.0*d+e)/(2.0*v) + (n+a)*f*(exp(e)-1.0);
        if(dL < 0.5){
            e_low = e;
        }else{
            e_high = e;
        }
        diff = fabs(dL-0.5);
    }
    return e*e;
}

void parse_argv(int argc,char** argv, string &in_file, string &gene_name_file, string &cell_name_file, string &in_file_extension, string &out_folder, int &N_threads, bool &print_extended_output, double &vmin, double &vmax, int &numbin, int &N_char, bool &no_norm, bool &max_v_output, bool &post_v_output){

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
			cout << VERSION << "\n";
			exit(0);
		}
	}

	int N_param(11);
	string extended_output("false");
	string no_norm_str("false");
	string max_v_str("true");
    string to_find[11][2] = {{"-f", "--file"},
							{"-d", "--destination"},
							{"-n", "--n_threads"},
							{"-e", "--extended_output"},
         					{"-vmin", "--variance_min"},
         					{"-vmax", "--variance_max"},
					    	{"-nbin", "--number_of_variance_bins"},
							{"-mtx_genes","--mtx_gene_name_file"},
							{"-mtx_cells","--mtx_cell_name_file"},
							{"-no_norm","--no_cell_size_normalization"},
							{"-max_v","--get_output_for_maxlik_variance"}};

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
				if(j==0) in_file = argv[idx+1];
				if(j==1) out_folder = argv[idx+1];
				if(j==2) N_threads = atoi(argv[idx+1]);
				if(j==3) extended_output = argv[idx+1];
				if(j==4) vmin = atof(argv[idx+1]);
				if(j==5) vmax = atof(argv[idx+1]);
				if(j==6) numbin = atoi(argv[idx+1]);
				if(j==7) gene_name_file = argv[idx+1];
				if(j==8) cell_name_file = argv[idx+1];
				if(j==9) no_norm_str = argv[idx+1];
				if(j==10) max_v_str = argv[idx+1];
				// add '/' to out_folder if not already
				if( j == 1 && out_folder.back() != '/' )
					out_folder = out_folder + '/';
            }
        }
		if (idx == 0 && j == 0){
        	cerr << "Error in argument parsing :\n"
            << "missing input file name\n";
			show_usage();
        }
    }

	if ( extended_output == "true" || extended_output == "1" )
		print_extended_output = true;

	if ( no_norm_str == "true" || no_norm_str == "1" )
		no_norm = true;
	
	if ( max_v_str == "true" || max_v_str == "1" || max_v_str == "only_max_output"){
		max_v_output = true;
		post_v_output = false;
		}
	// else if (max_v_str == "only_max_output"){
	// 	max_v_output = true;
	// 	post_v_output = false;
	// 	}


	// Get input file extension
	in_file_extension = (in_file.size() >= 3) ? in_file.substr(in_file.size() - 3) : in_file;
	for (char &ch : in_file_extension) {
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	cerr << "File type : " << in_file_extension << "\n";

	// Get number of Character in first row, for iobuffer
	string command = "head -1 " + in_file + "|wc -c>" + out_folder + "tmp";
	int out = system(command.c_str());
    ifstream myfile (out_folder + "tmp");
    string line;
    if (myfile.is_open()){
        getline (myfile,line);
        myfile.close();
		N_char = atoi(line.c_str());
		command = "rm " + out_folder + "tmp";
		out = system(command.c_str());
    }
    else{
		cerr << "Unable to open " + out_folder + "tmp\n";
		N_char = 10000000;
	}
	N_char = 2*N_char;
}

static void show_usage(void)
{
    cerr << "Usage: Sanity <option(s)> SOURCES\n"
         << "Options:\n"
         << "\t-h,--help\t\tShow this help message\n"
         << "\t-v,--version\t\tShow the current version\n"
         << "\t-f,--file\t\tSpecify the input transcript count text file (.mtx for Matrix Market File Format)\n"
		 << "\t-mtx_genes,--mtx_gene_name_file\tSpecify the gene name text file (only needed if .mtx input file)\n"
		 << "\t-mtx_cells,--mtx_cell_name_file\tSpecity the cell name text file (only needed if .mtx input file)\n"
         << "\t-d,--destination\tSpecify the destination path (default: pwd)\n"
         << "\t-n,--n_threads\t\tSpecify the number of threads to be used (default: 4)\n"
         << "\t-e,--extended_output\tOption to print extended output (default: false, choice: false,0,true,1)\n"
         << "\t-vmin,--variance_min\tMinimal value of variance in log transcription quotient (default: 0.001)\n"
         << "\t-vmax,--variance_max\tMaximal value of variance in log transcription quotient (default: 50)\n"
         << "\t-nbin,--number_of_bins\tNumber of bins for the variance in log transcription quotient  (default: 160)\n"
		 << "\t-no_norm,--no_cell_size_normalization\tOption to skip cell size normalization (default: false, choice: false,0,true,1)\n"
		 << "\t-max_v,--get_output_for_maxlik_variance\tOption to get the results for the prior variance (v_g) that maximizes the likelihood, \n\t"
		 <<	"\ti.e., without integrating over the posterior for v_g. (default: false, choice: false,0,true,1)\n";
	exit(0);
}

std::vector<double> fetch_row(int g, const std::string& in_file, const std::string& in_file_extension, const std::vector<RowBlock>& mtx_rows, const std::vector<std::streampos>& tsv_offsets, const int& C) {
	std::vector<double> row_data(C, 0.0);

	if (in_file_extension == "mtx") {
		// Fetch from Matrix Market format
		const RowBlock& row_block = mtx_rows[g];
		std::ifstream stream(in_file, std::ios::in | std::ios::binary);
        if (!stream) {
            throw std::runtime_error("RowReader: cannot open MTX " + in_file);
        }
		stream.seekg(row_block.offset);

		for (size_t idx = 0; idx < row_block.nnz; ++idx) {
			std::string line;
            std::getline(stream, line);
            std::istringstream iss(line);
            int g_idx, c_idx;
            double value;
            iss >> g_idx >> c_idx >> value;
            row_data[c_idx - 1] = value; // Convert to 0-based index
		}
	} else {
		// Fetch from TSV format
		std::ifstream infile(in_file, std::ios::in | std::ios::binary);
		if (!infile) {
            throw std::runtime_error("RowReader: cannot open TSV " + in_file);
        }

		infile.seekg(tsv_offsets[g]);
		std::string line;
		std::getline(infile, line);
		std::istringstream iss(line);
		std::string token;
		std::string gene_id;
		std::getline(iss, gene_id, '\t');  // skip gene name
		for (int c = 0; c < C && std::getline(iss, token, '\t'); ++c) {
    		row_data[c] = std::stod(token);
		}
		infile.close();
	}

	return row_data;
}
