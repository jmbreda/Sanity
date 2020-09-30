clear all;

%%%%% Get the raw dataset on Zenedo: https://doi.org/10.5281/zenodo.3996271 %%%%

Datasets = {'Gruen_ESC_SC_2i','Gruen_ESC_SC_serum','Gruen_ESC_RNA_2i','Gruen_ESC_RNA_serum','Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo' 'Simulated_Baron_Independent_Genes','Simulated_Branched_Random_Walk'};

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};

% Get linear expression from different normalizations
for d = Datasets
	for my_norm = setdiff(My_norm,'sctransform');

		% Read input table
		if strcmp(my_norm{:},'Deconvolution')
			in_file = ['data/' d{:} '_UMI_counts.txt'];
		else
			in_file = ['data/' d{:} '_' my_norm{:} '_normalization.txt'];
		end
		T = readtable(in_file,'ReadRowNames',1,'delimiter','\t');

		% Output files
		out_file = ['data/' d{:} '_' my_norm{:} '_normalization_lin.mat'];

		% Get linear expression from different normalizations
		switch my_norm{:}
		case 'RawCounts'
			M = exp(T{:,:})-1;
		case 'TPM'
			M = exp(T{:,:})-1;
		case 'DCA'
			M = T{:,:};
		case 'Deconvolution'
			sizeFactor = readtable(['data/' d{:} '_' my_norm{:} '_normalization.txt'])';
			M = T{:,:}./sizeFactor{:,:};
		case 'MAGIC'
			M = T{:,:};
		case 'Sanity'
			M = exp(T{:,:});
		case 'SAVER'
			M = T{:,:};
		case 'scImpute'
			M = T{:,:};
		case 'scVI'
			M = T{:,:};
		end

		% Save matrix
		save(out_file,'M','-v7.3');
	end
end

% Get log expression from different normalizations
for d = Datasets
	for my_norm = My_norm

		% Read input table
		if strcmp(d{:},'Deconvolution')
			in_file = ['data/' d{:} '_UMI_counts.txt'];
		else
			in_file = ['data/' d{:} '_' my_norm{:} '_normalization.txt'];
		end
		T = readtable(in_file,'ReadRowNames',1,'delimiter','\t');

		% Output files
		out_file = ['data/' d{:} '_' my_norm{:} '_normalization.mat'];

		% Get log expression from different normalizations
		switch my_norm{:}
		case 'RawCounts'
			M = T{:,:};
		case 'TPM'
			M = T{:,:};
		case 'DCA'
			M = log(T{:,:});
		case 'Deconvolution'
			sizeFactor = readtable(['data/' d{:} '_' my_norm{:} '_normalization.txt'])';
			M = log(T{:,:}./sizeFactor{:,:}+1);
		case 'MAGIC'
			M = T{:,:};
			M = max(0,M);
			if min(min(M)) == 0
				p=1;
			else
				p=0;
			end
			M = log(M+p);
		case 'Sanity'
			M = T{:,:};
		case 'SAVER'
			M = log(T{:,:});
		case 'scImpute'
			M = log(T{:,:}+1);
		case 'sctransform'
			M = T{:,:};
		case 'scVI'
			M = log(T{:,:});
		end

		% Save matrices
		save(out_file,'M','-v7.3');
	end
end

