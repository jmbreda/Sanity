#include <ReadInputFiles.h>

void Get_G_C_UMIcountMatrix(string in_file, int &N_rows, int &G, int &C, int N_char){

    int g, c;
    char *ss = new char [N_char];
    char *sc = new char [N_char];
	
	// Open input file
    FILE *infp;
    infp = (FILE *) fopen(in_file.c_str(),"r");
    if(infp == NULL){
        fprintf(stderr,"Cannot open input file %s\n",in_file.c_str());
        exit(EXIT_FAILURE);
    }

    // Count cell. First line should have the names of the columns (cell names)
    C = 0;
    fgets(ss,N_char,infp);
    strcpy(sc,ss);
    char *token = strtok(ss," \t");
    while(token){
        ++C;
        token = strtok(NULL," \t");
    }
    C = C-1;

	// Count rows
    N_rows = 0;
    while( fgets(ss,N_char,infp) ){
        ++N_rows;
    }
    fclose(infp);

    // Total Count per gene :
    double n(0);
    G = 0;

    // Reopen file to read the counts per gene
    infp = (FILE *) fopen(in_file.c_str(),"r");
    /***first line has names***/
    fgets(ss,N_char,infp);
    char *retval;
    for(g=0; g<N_rows; ++g){
        retval = fgets(ss,N_char,infp);
        if(retval == NULL){
            fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
            exit(EXIT_FAILURE);
        }
        /***cut the line based on spaces/tabs***/
        strcpy(sc,ss);
        token = strtok(ss," \t");
        n = 0;
        for(c=0; c<C; ++c){
            /****go to the next in the list of words****/
            token = strtok(NULL," \t");
            if(token != NULL){
                n += atoi(token);/**total count this gene**/
            }else{
                fprintf(stderr,"Error: not enough fields on line number %d:\n%s\n",g,sc);
                exit(EXIT_FAILURE);
            }
        }
		// Add gene if total count bigger than 0
        if (n > 0){
            G++;
        }
        token = strtok(NULL," \t");
        if(token != NULL){
            fprintf(stderr,"Error: too many fields on line number %d:\n%s\n",g,sc);
        }
    }

    fclose(infp);
    delete[] ss;
    delete[] sc;

	return;
}

void ReadUMIcountMatrix(string in_file, double **n_c, double *N_c, double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, int N_char){

    int row, g, c;
    char *ss = new char [N_char];
    char *sc = new char [N_char];

    FILE *infp;
    infp = (FILE *) fopen(in_file.c_str(),"r");
    if(infp == NULL){
        fprintf(stderr,"Cannot open input file %s\n",in_file.c_str());
        exit(EXIT_FAILURE);
    }

    /***First line should have the names of the columns (sample names)***/
    fgets(ss,N_char,infp);

    /***now go fill in the cell names***/
	char *token = strtok(ss," \t");
	// Ignore first word
    token = strtok(NULL," \t");
    c=0;
    int tmp;
    while(token){
        cell_names[c] = string(token);
        tmp = cell_names[c].find('\n');
        if (tmp > -1){
            cell_names[c].erase(tmp,1);
        }
        ++c;
        token = strtok(NULL," \t");
    }

	
    // Temporary count per gene and cell
    double *n_c_tmp = new double [C];
    // Temporary Count per gene :
	double n_tmp;
	// temporary gene name
	string gene_names_tmp;
	// expressed gene index
	g = 0;

    /**reopen file to read the counts***/
    infp = (FILE *) fopen(in_file.c_str(),"r");
    // First line has names
    fgets(ss,N_char,infp);
    char *retval;
    for(row=0; row<N_rows; ++row){
        retval = fgets(ss,N_char,infp);
        if(retval == NULL){
            fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
            exit(EXIT_FAILURE);
        }
        /***cut the line based on spaces/tabs***/
        strcpy(sc,ss);
        token = strtok(ss," \t");

        gene_names_tmp = token;
        n_tmp = 0;
        for(c=0; c<C; ++c){
            /****go to the next in the list of words****/
            token = strtok(NULL," \t");
            if(token != NULL){
                n_c_tmp[c] = atoi(token);
                n_tmp += n_c_tmp[c];// total count this gene
                N_c[c] += n_c_tmp[c];// total count for this cell
            }else{
                fprintf(stderr,"Error: not enough fields on line number %d:\n%s\n",g,sc);
                exit(EXIT_FAILURE);
            }
        }
		// if total count for the gene not 0, copy gene_name, gene count and total count
        if (n_tmp > 0){
			n[g] = n_tmp;
			gene_names[g] = gene_names_tmp;
			for(c=0; c<C; ++c)
				n_c[g][c] = n_c_tmp[c];
			g++;
        }
        token = strtok(NULL," \t");
        if(token != NULL){
            fprintf(stderr,"Error: too many fields on line number %d:\n%s\n",g,sc);
        }
    }

    fclose(infp);
    delete[] ss;
    delete[] sc;
    delete[] n_c_tmp;
	
	return;
}

void Get_G_C_MTX(string in_file, int &N_rows, int &G, int &C, map<int,int> &gene_idx){

	// declare variables
    FILE * infp;
	char ss[1024];
	char * retval = NULL;
	char * token = NULL;

	// Read mtx file
    infp = fopen(in_file.c_str(), "r");
    if (infp == NULL){
        fprintf(stderr,"Cannot open input file %s\n",in_file.c_str());
        exit(EXIT_FAILURE);
	}

	// Ignore header and read first row (G C totUMI)
	while ( (retval = fgets(ss,1024,infp)) != NULL ){
		if ( retval[0] == '%' ){
			continue;
		}else{
			token = strtok(retval," ");
			// Highest number of genes expressed and not
			N_rows = atoi(token);
			token = strtok(NULL," \t");
			// Number of cells
			C = atoi(token);
			break;
		}
	}

	// Boolean vector of expressed genes, start with all gene not expressed.
	bool expressed_genes[N_rows];
	int g;
	for (g=0;g<N_rows;g++)
		expressed_genes[g] = 0;
	
	// Read following rows rest of
    while ( (retval = fgets(ss,1024,infp) ) != NULL) {
		// Read values gene idx and add as expressed
		token = strtok(retval," \t");
		g = atoi(token);
		g=g-1;
		// Update gene_idx map and count gene if not seen before
		if( ! expressed_genes[g] ){
			expressed_genes[g] = 1;
			G++;
        }
    }
    fclose(infp);

	// Create gene_idx map
	int idx = 0;
	for(g=0;g<N_rows;g++){
		if(expressed_genes[g]){
			gene_idx[g] = idx++;
		}else{
			gene_idx[g] = -1;
		}
	}

	return;
}

void ReadMTX(string mtx_file, string gene_name_file, string cell_name_file, double **n_c, double *N_c, double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, map<int,int> gene_idx){

	// declare variables
    FILE * infp;
	char ss[1024];
	char * retval = NULL;
	char * token = NULL;
	int g,c,count,tmp;

	// Read gene_names
	if ( gene_name_file != "none" ){
		infp = fopen(gene_name_file.c_str(), "r");
		if (infp == NULL){
			fprintf(stderr,"Cannot open input file %s\n",gene_name_file.c_str());
			exit(EXIT_FAILURE);
		}

		g = 0;	
		while ( (retval = fgets(ss,1024,infp)) != NULL ){
			if ( retval[0] != '%' ){
				if (gene_idx[g] != -1){
					token = strtok(retval," \t");
					gene_names[gene_idx[g]] = string(token);
					token = strtok(NULL," \t");
					if (token != NULL )
						gene_names[gene_idx[g]] = gene_names[gene_idx[g]] + "|" + string(token);

					tmp = gene_names[gene_idx[g]].find('\n');
					if (tmp > -1)
						gene_names[gene_idx[g]].erase(tmp,1);
				}
				++g;
			}
		}
	}else{
		fprintf(stdout,"No gene name file\n");
		for(g=0;g<G;g++)
			gene_names[g] = "Gene_" + to_string(g);
	}

	// Read cell names
	if ( cell_name_file != "none" ){
		infp = fopen(cell_name_file.c_str(), "r");
		if (infp == NULL){
			fprintf(stderr,"Cannot open input file %s\n",cell_name_file.c_str());
			exit(EXIT_FAILURE);
		}

		c = 0;
		while ( (retval = fgets(ss,1024,infp)) != NULL ){
			if ( retval[0] != '%' ){
				token = strtok(retval," \t");
				cell_names[c] = string(token);
				
				tmp = cell_names[c].find('\n');
				if (tmp > -1){
					cell_names[c].erase(tmp,1);
				}
				++c;
			}
		}
	}else{
		fprintf(stdout,"No cell name file\n");
		for(c=0;c<C;c++)
			cell_names[c] = "Cell_" + to_string(c);
	}


	// Read mtx file
    infp = fopen(mtx_file.c_str(), "r");
    if (infp == NULL){
		fprintf(stderr,"Cannot open input file %s\n",mtx_file.c_str());
        exit(EXIT_FAILURE);
	}

	// Ignore header and first row (G C totUMI)
	while ( (retval = fgets(ss,1024,infp)) != NULL ){
		if ( retval[0] == '%' ){
			continue;
		}else{
			token = strtok(retval," ");
			break;
		}
	}

	// Read following rows rest of
    while ( (retval = fgets(ss,1024,infp) ) != NULL) {
		// Read values gene idx and add as expressed
		token = strtok(retval," ");
		g = atoi(token);
		token = strtok(NULL," ");
		c = atoi(token);
		token = strtok(NULL," ");
		count = atoi(token);

		g = g-1;
		c = c-1;
		// fill in the count matrix
		n_c[gene_idx[g]][c] = count;
		n[gene_idx[g]] += count;
		N_c[c] += count;
    }
    fclose(infp);

	return;
}
