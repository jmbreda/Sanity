#include <ReadInputFiles.h>

void Get_G_C_UMIcountMatrix(string in_file,
							int &N_rows,
							int &G,
							int &C,
							int N_char,
							std::vector<std::streampos> &tsv_offsets,
							std::vector<double> &N_c,
							std::vector<double> &n,
							std::vector<std::string> &cell_names,
							std::vector<std::string> &gene_names){

	std::ifstream infp(in_file.c_str(), std::ios::in | std::ios::binary);
	if(!infp.is_open()){
		fprintf(stderr,"Cannot open input file %s\n",in_file.c_str());
		exit(EXIT_FAILURE);
	}

	std::vector<char> ss(static_cast<std::size_t>(N_char));
	std::vector<char> sc(static_cast<std::size_t>(N_char));

	// Count cell. First line should have the names of the columns (cell names)
	if(!infp.getline(ss.data(), N_char)){
		fprintf(stderr,"Error: Unable to read header line from %s\n", in_file.c_str());
		exit(EXIT_FAILURE);
	}
	std::strncpy(sc.data(), ss.data(), static_cast<std::size_t>(N_char));
	sc[N_char-1] = '\0';
	char *token = strtok(ss.data()," \t,");
	cell_names.clear();
	if (token == NULL){
		fprintf(stderr,"Error: header line missing first field in %s\n", in_file.c_str());
		exit(EXIT_FAILURE);
	}
	// remaining tokens are cell names
	token = strtok(NULL," \t,");
	while(token){
		std::string name(token);
		size_t pos = name.find('\r');
		if (pos != std::string::npos){
			name.erase(pos,1);
		}
		pos = name.find('\n');
		if (pos != std::string::npos){
			name.erase(pos,1);
		}
		cell_names.push_back(name);
		token = strtok(NULL," \t,");
	}
	C = static_cast<int>(cell_names.size());
	std::vector<double> cell_totals(static_cast<std::size_t>(C), 0.0);
	std::vector<double> gene_totals;
	gene_names.clear();

	// Prepare to record offsets and counts
	N_rows = 0;
	G = 0;
	tsv_offsets.clear();

	while(true){
		std::streampos row_offset = infp.tellg();
		if(!infp.getline(ss.data(), N_char)){
			break;
		}
		++N_rows;
		std::strncpy(sc.data(), ss.data(), static_cast<std::size_t>(N_char));
		sc[N_char-1] = '\0';

		token = strtok(ss.data()," \t,");
		if(token == NULL){
			continue;
		}
		std::string gene_id(token);
		size_t pos = gene_id.find('\r');
		if (pos != std::string::npos){
			gene_id.erase(pos,1);
		}
		pos = gene_id.find('\n');
		if (pos != std::string::npos){
			gene_id.erase(pos,1);
		}

		double row_sum = 0.0;
		for(int c=0; c<C; ++c){
			token = strtok(NULL," \t,");
			if(token != NULL){
				double value = stod(token);/**total count this gene**/
				row_sum += value;
				cell_totals[static_cast<std::size_t>(c)] += value;
			}else{
				fprintf(stderr,"Error: not enough fields on line number %d:\n%s\n",N_rows-1,sc.data());
				exit(EXIT_FAILURE);
			}
		}
		// Add gene if total count bigger than 0
		if (row_sum > 0){
			G++;
			tsv_offsets.push_back(row_offset);
			gene_totals.push_back(row_sum);
			gene_names.push_back(gene_id);
		}
		token = strtok(NULL," \t,");
		if(token != NULL){
			fprintf(stderr,"Error: too many fields on line number %d:\n%s\n",N_rows-1,sc.data());
		}
	}

	infp.close();

	N_c.assign(cell_totals.begin(), cell_totals.end());
	n.assign(gene_totals.begin(), gene_totals.end());

	return;
}

void ReadUMIcountMatrix(string in_file, double **n_c, double *N_c, double *n, string *gene_names, string *cell_names, int N_rows, int G, int C, int N_char){

    int row, g, c;
    char *retval;
    char *ss = new char [N_char];
    char *sc = new char [N_char];

    FILE *infp;
    infp = (FILE *) fopen(in_file.c_str(),"r");
    if(infp == NULL){
        fprintf(stderr,"Cannot open input file %s\n",in_file.c_str());
        exit(EXIT_FAILURE);
    }

    /***First line should have the names of the columns (sample names)***/
    retval = fgets(ss,N_char,infp);

    /***now go fill in the cell names***/
	char *token = strtok(ss," \t,");
	// Ignore first word
    token = strtok(NULL," \t,");
    c=0;
    int tmp;
    while(token){
        cell_names[c] = string(token);
        tmp = cell_names[c].find('\n');
        if (tmp > -1){
            cell_names[c].erase(tmp,1);
        }
        ++c;
        token = strtok(NULL," \t,");
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
    retval = fgets(ss,N_char,infp);
    for(row=0; row<N_rows; ++row){
        retval = fgets(ss,N_char,infp);
        if(retval == NULL){
            fprintf(stderr,"Error: Couldn't read a line at row number %d\n",g);
            exit(EXIT_FAILURE);
        }
        /***cut the line based on spaces/tabs***/
        strcpy(sc,ss);
        token = strtok(ss," \t,");

        gene_names_tmp = token;
        n_tmp = 0;
        for(c=0; c<C; ++c){
            /****go to the next in the list of words****/
            token = strtok(NULL," \t,");
            if(token != NULL){
                n_c_tmp[c] = stod(token);
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
        token = strtok(NULL," \t,");
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

void Get_G_C_MTX(string in_file, int &N_rows, int &G, int &C, map<int,int> &gene_idx, vector<RowBlock> &mtx_rows, vector<double> &N_c, vector<double> &n){

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
			N_rows = stoi(token);
			token = strtok(NULL," \t");
			// Number of cells
			C = stoi(token);
			break;
		}
	}

	// Boolean vector of expressed genes, start with all gene not expressed.
	vector<char> expressed_genes(N_rows, 0);
	vector<RowBlock> row_blocks_tmp(N_rows);
	vector<double> gene_totals(N_rows, 0.0);
	vector<double> cell_totals(static_cast<std::size_t>(C), 0.0);
	int g;
	for (g=0;g<N_rows;g++) {
		row_blocks_tmp[g].offset = -1;
		row_blocks_tmp[g].nnz = 0;
	}
	
	// Read following rows rest of
	G=0;
	while ( true ) {
		long offset = ftell(infp);
		retval = fgets(ss,1024,infp);
		if ( retval == NULL ){
			break;
		}
		if ( retval[0] == '%' || retval[0] == '\n' || retval[0] == '\0'){
			continue;
		}
		// Read values gene idx and add as expressed
		token = strtok(retval," \t\r\n");
		if (token == NULL){
			continue;
		}
		g = stoi(token);
		g=g-1;
		if ( g < 0 || g >= N_rows ){
			fprintf(stderr,"Error: gene index %d out of range in MTX file %s\n", g+1, in_file.c_str());
			fclose(infp);
			exit(EXIT_FAILURE);
		}
		token = strtok(NULL," \t\r\n");
		if (token == NULL){
			fprintf(stderr,"Error: missing cell index for gene %d in MTX file %s\n", g+1, in_file.c_str());
			fclose(infp);
			exit(EXIT_FAILURE);
		}
		int c = stoi(token) - 1;
		if ( c < 0 || c >= C ){
			fprintf(stderr,"Error: cell index %d out of range in MTX file %s\n", c+1, in_file.c_str());
			fclose(infp);
			exit(EXIT_FAILURE);
		}
		token = strtok(NULL," \t\r\n");
		if (token == NULL){
			fprintf(stderr,"Error: missing count for gene %d cell %d in MTX file %s\n", g+1, c+1, in_file.c_str());
			fclose(infp);
			exit(EXIT_FAILURE);
		}
		double count = stod(token);
		// Update gene_idx map and count gene if not seen before
		if( ! expressed_genes[g] ){
			expressed_genes[g] = 1;
			G++;
			row_blocks_tmp[g].offset = static_cast<long long>(offset);
			row_blocks_tmp[g].nnz = 1;
		}else{
			row_blocks_tmp[g].nnz += 1;
        }
		gene_totals[g] += count;
		cell_totals[static_cast<std::size_t>(c)] += count;
    }
    fclose(infp);

	// Create gene_idx map
	int idx = 0;
	mtx_rows.clear();
	mtx_rows.resize(G);
	n.assign(static_cast<std::size_t>(G), 0.0);
	for(g=0;g<N_rows;g++){
		if(expressed_genes[g]){
			gene_idx[g] = idx;
			mtx_rows[idx] = row_blocks_tmp[g];
			n[static_cast<std::size_t>(idx)] = gene_totals[g];
			idx++;
		}else{
			gene_idx[g] = -1;
		}
	}
	N_c.assign(cell_totals.begin(), cell_totals.end());

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
		g = stoi(token);
		token = strtok(NULL," ");
		c = stoi(token);
		token = strtok(NULL," ");
		count = stod(token);

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

std::vector<std::string> Read_CellNames(const std::string &filename){
	std::vector<std::string> cell_names;
	if (filename == "none" || filename.empty()){
		return cell_names;
	}
	std::ifstream input(filename.c_str(), std::ios::in);
	if (!input.is_open()){
		fprintf(stderr,"Cannot open input file %s\n", filename.c_str());
		exit(EXIT_FAILURE);
	}
	std::string line;
	while (std::getline(input, line)){
		if (line.empty() || line[0] == '%'){
			continue;
		}
		if (!line.empty() && line.back() == '\r'){
			line.pop_back();
		}
		cell_names.push_back(line);
	}
	return cell_names;
}

std::vector<std::string> Read_GeneNames(const std::string &filename){
	std::vector<std::string> gene_names;
	if (filename == "none" || filename.empty()){
		return gene_names;
	}
	std::ifstream input(filename.c_str(), std::ios::in);
	if (!input.is_open()){
		fprintf(stderr,"Cannot open input file %s\n", filename.c_str());
		exit(EXIT_FAILURE);
	}
	std::string line;
	while (std::getline(input, line)){
		if (line.empty() || line[0] == '%'){
			continue;
		}
		if (!line.empty() && line.back() == '\r'){
			line.pop_back();
		}
		gene_names.push_back(line);
	}
	return gene_names;
}
