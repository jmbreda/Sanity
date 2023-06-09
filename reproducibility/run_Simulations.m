clc; clear all;

sim_seed = 42;
rng(sim_seed);

% load Baron raw UMI counts
T = readtable('data/Baron_UMI_counts.txt','ReadRowNames',1,'delimiter','\t');

% Loop over the 2 simulations
for sim_name = {'Simulated_Baron_Independent_Genes','Simulated_Branched_Random_Walk'} 

	% get dimensions
	[N_gene,N_cell] = size(T{:,:});

	% Get total UMI per cell, randomly shuffle
	N_c = sum(T{:,:},1);
	N_c = N_c(randperm(N_cell));

	% Get gene mean transcription quotient from raw count Baron
	mu_tilde_g = log(sum(T{:,:},2)./sum(sum(T{:,:})));

	% Simulate gene variance
	sig2_g = exprnd(2,N_gene,1);

	switch sim_name{:}
	case 'Simulated_Baron_Independent_Genes'

		%1. get delta Norm(0,1)
		delta_true = normrnd(0,1,N_gene,N_cell);

		% 2. Shift <delta_g> = 0 
		%    Multiply by lamdba_g s.t. var(lambda_g*delta_g) = sig2_g => lambda = sqrt( sig2_g/var(delta_g));
		lambda = sqrt( sig2_g./var(delta_true,0,2) );
		delta_true = lambda.*(delta_true-mean(delta_true,2));

		% 3. add mean s.t sum_g exp( mu_g + delta_gc ) = 1
		mu_g = mu_tilde_g - sig2_g/2;
		E = mu_g+delta_true;

		neighbor = 0:N_cell-1;

	case 'Simulated_Branched_Random_Walk'

		N_path = 149;
		length_path = 13;

		% 1. get delta Norm(0,1)
		c = 1;
		neighbor = [];
		delta_true = [];
		for k = 1:N_path
			if k==1
				neighbor(c) = 0;
				d_0 = zeros(N_gene,1);
			else
				neighbor(c) = randi(c);
				d_0 = delta_true(:,neighbor(c));
			end
			for i = 1:length_path
				if i==1
					delta_true(:,c) = d_0 + normrnd(0,1,N_gene,1);
				else
					neighbor(c) = c-1;
					delta_true(:,c) = delta_true(:,c-1) + normrnd(0,1,N_gene,1);
				end
				c = c+1;
			end
		end

		% 2. Shift <delta_g> = 0 
		%    And multiply by lamdba_g s.t. var(lambda_g*delta_g) = sig2_g => lambda = sqrt( sig2_g/var(delta_g));
		lambda = sqrt( sig2_g./var(delta_true,0,2) );
		delta_true = lambda.*(delta_true-mean(delta_true,2));

		% 3. add mean s.t sum_g exp( mu_g + delta_gc ) = 1
		mu_g = mu_tilde_g - sig2_g/2;
		E = mu_g+delta_true;
	end

	% Take poisson sampling
    UMI = poissrnd( N_c.*exp(E) );

	% Remove non expressed genes
    idx_0 = find(sum(UMI,2)==0);
    if ~isempty(idx_0)
        mu_g(idx_0) = [];
        sig2_g(idx_0) = [];
        delta_true(idx_0,:) = [];
        E(idx_0,:) = [];
        UMI(idx_0,:) = [];
        N_gene = size(mu_g,1);
    end

	% Compute cell-cell euclidean distances
    D_true = pdist(E')';

	% Make UMI count table
    Gene_name = [];
    Cell_name = [];
    for i=1:N_gene
        Gene_name{i,1} = ['Gene_' num2str(i)];
    end
    for i=1:N_cell
        Cell_name{i,1} = ['Cell_' num2str(i)];
    end
    UMI_table = array2table(UMI);
    UMI_table.Properties.RowNames = Gene_name;
    UMI_table.Properties.VariableNames = Cell_name;

	% Save the table and workspace
    writetable(UMI_table,['data/' sim_name{:} '_UMI_counts.txt'],'writeRowNames',1,'delimiter','\t');
	save(['data/' sim_name{:} '.mat'])

    %dlmwrite(['data/distance_true.txt'],D_true,'delimiter','\t');
    %dlmwrite(['data/delta_true.txt'],delta_true,'delimiter','\t');
    %dlmwrite(['data/ltq_true.txt'],E,'delimiter','\t');
    %dlmwrite(['/variance_true.txt'],sig2_g,'delimiter','\t');
end
