close all; clear all;

% Initalize random number generator
sim_seed = 42;
rng(sim_seed);

% Number of cells
N_cell = 1000;

% 0. Simulate 3 lognormally distributed genes
mu = [-1.9 -2 -2.8]';
sig2 = [0 .1 1]';
e = normrnd(repmat(mu,1,N_cell),repmat(sqrt(sig2),1,N_cell),3,N_cell);

% Take Poisson sampling form the expression values
T_capt = poissrnd(exp(e));

% Simulate total count per cell
N_c = round(10.^normrnd(4,.2,1,N_cell));

% Add a gene to get to the desired total gene per cell
T_capt(4,:) = N_c - sum(T_capt,1);

% Make and print gene UMI count table
dataset = ['Low_expression'];
out_dir = ['output/' dataset];
mkdir(out_dir);
N_gene = size(T_capt,1);
for i=1:N_gene
    Gene{i,1} = ['Gene_' num2str(i)];
end
Transcript_captured = [cell2table(Gene) array2table(T_capt)];
writetable(Transcript_captured,[out_dir '/Transcript_captured.txt'],'delimiter','\t')

% Run Sanity
my_file = [out_dir '/Transcript_captured.txt'];
[stat,out] = system(['Sanity -f ' my_file ' -d ' out_dir ' -e 1 -n 4 -vmin 0.001 -vmax 50 -nbin 160'])

% Load gene variance and variance likelihood
v = load([out_dir '/variance.txt']);
L = readtable([out_dir '/likelihood.txt'],'ReadRowNames',1,'ReadVariableName',0);

% variance bins in first row
v_bin = L{1,:};

% Likelihood in subsequent rows
L = L{2:4,:};

% Load Sanity outputs
ltq = readtable([out_dir '/log_transcription_quotients.txt'],'ReadRowNames',1);
d_ltq = readtable([out_dir '/ltq_error_bars.txt'],'ReadRowNames',1);
ltq = ltq{:,:};
d_ltq = d_ltq{:,:};

% Take only the 3 lognormally distributed genes
T_capt = T_capt(1:3,:);

% Define colors
my_colors = lines(3);

% Plot
figure('visible','off'); f=1;
for g = 1:3
	% Distribution of transcription activity ofa each gene
	subplot(3,3,f);f=f+1;
	[h,x] = hist(e(g,:),30);
	h = bar(x,h,1);
	set(h,'FaceColor',my_colors(g,:),'linestyle','none')
	xlim([-6 2])
	xlabel('True transcription activity')
	ylabel('Frequency')
	title(['\mu = ' num2str(mu(g)) ' \sigma^2=' num2str(sig2(g))],'FontWeight','normal')
end

for g = 1:3
	% Inferred likelihood of variance of each gene
	subplot(3,3,f);f=f+1;
	[h,x] = hist(e(g,:),30);
	semilogx(v_bin,L(g,:),'color',my_colors(g,:),'linewidth',2);
    hold on
    plot([v(g) v(g)],[0 max(L(g,:))],'k-')
	xlabel('variance')
	ylabel('Likelihood')
    xlim([1e-3 50])
    set(gca,'xtick',10.^[-3:1])
end

% Histogram of Simulated UMI count for each gene
subplot(3,3,f);f=f+1;
bins = 0:max(T_capt(:));
h = histc(T_capt',bins);
bar(bins,h,1)
set(gca,'yscale','log','ytick',10.^[0:3],'yticklabel',{'1' '10' '100' '1000'})
xlim([-0.5,max(bins)+.5])
ylim([0.5 2000])
xlabel('UMI count')

% Print
dim = [18 15];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S21','-dpdf');
