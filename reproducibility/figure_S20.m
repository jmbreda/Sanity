close all; clear all;

% Initalize random number generator
sim_seed = 42;
rng(sim_seed);

% Number of cells
N_cell = 2000;

% 0. Simulate lognormal distributed gene
mu = log(1);
sig2 = 1;
e(1,:) = normrnd(mu,sqrt(sig2),1,N_cell);

% 1. Simulate Bimodal distributed gene
mu = [log(1) log(20)];
sig2 = [.1 .1];
e(2,:) = [ normrnd(mu(1),sqrt(sig2(1)),1,N_cell/2) normrnd(mu(2),sqrt(sig2(2)),1,N_cell/2)];

% 2. Simulate Uniformly distributed gene
e(3,:) = unifrnd(log(1),log(exp(4)),1,N_cell);

% Take Poisson sampling form the expression values
T_capt = poissrnd(exp(e));

% Simulate total count per cell
N_c = round(10.^normrnd(4,.2,1,N_cell));

% 3. Add a gene to get to the desired total gene per cell
T_capt(4,:) = N_c - sum(T_capt,1);

% Get true log transcription quotients
e = e - log(N_c);

% Make and print gene UMI count table
dataset = ['Norm_Bimodal_Unif'];
out_dir = ['output/' dataset];
mkdir(out_dir);
N_gene = size(T_capt,1);
for i=1:N_gene
	Gene{i,1} = ['Gene_' num2str(i)];
end
Transcript_captured = [cell2table(Gene) array2table(T_capt)];
writetable(Transcript_captured,[out_dir '/Transcript_captured.txt'],'delimiter','\t')

% Run Sanity on the simulated UMI count table
my_file = [out_dir '/Transcript_captured.txt'];
[stat,out] = system(['Sanity -f ' my_file ' -d ' out_dir ' -e 1 -n 4'])

% Load Sanity outputs
ltq = readtable([out_dir '/log_transcription_quotients.txt'],'ReadRowNames',1);
d_ltq = readtable([out_dir '/ltq_error_bars.txt'],'ReadRowNames',1);
ltq = ltq{:,:};
d_ltq = d_ltq{:,:};

% Compute TPM
tpm = log( T_capt./sum(T_capt,1)*median(sum(T_capt,1)) + 1);

% Plot
figure('visible','off'); f=1;
for g = 1:3

	% True ltq distribution
	subplot(3,5,f);f=f+1;
	[h,x] = hist(e(g,:),30);
	bar(x,h,1,'k')
	if g==3
		xlabel('True ltq')
	end
	ylabel('Frequency')

	% UMI count distribution
	subplot(3,5,f);f=f+1;
	bins = 0:max(T_capt(g,:));
	h = histc(T_capt(g,:),bins);
	bar(bins,h,1,'k')
	if g==3
		xlabel('UMI count')
	end

	% log TPM distribution
	subplot(3,5,f);f=f+1;
	[h,x] = hist(tpm(g,:),30);
	bar(x,h,1,'k')
	if g==3
		xlabel('log TPM')
	end

	% Sanity inferred ltq distribution
	subplot(3,5,f);f=f+1;
	[h,x] = hist(ltq(g,:),30);
	bar(x,h,1,'k')
	if g==3
		xlabel('Sanity ltq')
	end

	% Sanity posterior distribution by adding small gaussian contribution of each cell
	subplot(3,5,f);f=f+1;
	n_bin = 100;
	x = linspace(min(ltq(g,:)-2*d_ltq(g,:)),max(ltq(g,:)+2*d_ltq(g,:)),n_bin);
	dx = x(2)-x(1);
	p_x = zeros(1,n_bin);
	for c = 1:N_cell
		p_x = p_x + normpdf(x,ltq(g,c),d_ltq(g,c))*dx;
	end
	p_x = p_x/N_cell;
	plot(x,p_x,'k')
	if g==3
		xlabel('Sanity ltq posterior prob.')
	end
end

% Print
dim = [30 12];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S20','-dpdf');
