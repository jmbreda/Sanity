clc; clear all;

sim_seed = 42;
rng(sim_seed);

T = readtable('data/Baron_UMI_counts.txt','ReadRowNames',1,'delimiter','\t');
[N_gene,N_cell] = size(T{:,:});

% Total UMI per cell
N_c = sum(T{:,:},1);

% mean and variance in expressions levels from Baron
load('data/Baron_Sanity_normalization.mat');
mu_g = mean(M,2);
sig2_g = unifrnd(0,6,N_gene,1);

% Simulate lognormal expression
E = median(N_c)*exp(normrnd(repmat(mu_g,1,N_cell),repmat(sqrt(sig2_g),1,N_cell)));
% Add fluctuations in total UMI count per cell
e = bsxfun(@times,bsxfun(@rdivide,E,sum(E,1)),N_c);
% Add Poisson noise
T_capt = poissrnd(e);

% Remove none expressed genes
ind_0 = find(sum(T_capt,2)==0);
if ~isempty(ind_0) 
    E(ind_0,:) = [];
    e(ind_0,:) = [];
    N_gene = size(E,1);
    T_capt(ind_0,:) = [];
end

for i=1:N_gene
	Gene{i,1} = ['Gene_' num2str(i)];
end
Transcript_captured = [cell2table(Gene) array2table(T_capt)];

out_dir = ['data/SimulatedBaron'];
if ~exist(out_dir,'dir')
	mkdir(out_dir)
end

% Save the workspace in my_sim.mat
save([out_dir '/my_sim'])

% Print simulated umi count matrix
writetable(Transcript_captured,[out_dir '/UMI_counts.txt'],'delimiter','\t')
