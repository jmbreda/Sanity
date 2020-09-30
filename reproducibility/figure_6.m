clear all; close all; clc;

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};

load('data/my_colors.txt');
my_colors = my_colors(1:10,:);
% Set default colors to my_colors
set(gca,'defaultAxesColorOrder',my_colors);

% Get true mean, variance and ltq (file from run_Simulations.m)
load(['data/Simulated_Baron_Independent_Genes.mat']);
%True.mean_tq = exp(mu_g)/sum(exp(mu_g));
True.mean_tq = mean(exp(E),2)/sum(mean(exp(E),2),1);
%True.mean_tq = mean(E./sum(E,2),2);
True.var = sig2_g;
%True.var = var(log(E),0,2);
True.lfc = E - mean(E,2);

% Get mean, variance and ltq from different normalisations
FQ = [];
for my_norm = My_norm
	switch my_norm{:}
	case 'sctransform'
		mu = load(['data/Simulated_Baron_Independent_Genes_sctransform_mean_expression.txt']);
		FQ.mean_tq.(my_norm{:}) = mu/sum(mu);
	otherwise
		load(['data/Simulated_Baron_Independent_Genes_' my_norm{:} '_normalization_lin.mat']);
		FQ.mean_tq.(my_norm{:}) = mean(M,2)/sum(mean(M,2),1);
	end
	load(['data/Simulated_Baron_Independent_Genes_' my_norm{:} '_normalization.mat']);
	FQ.var.(my_norm{:}) = var(M,0,2);
	FQ.lfc.(my_norm{:}) = M - mean(M,2);
end
T = readtable(['data/Simulated_Baron_Independent_Genes_Sanity_variance.txt']);
FQ.var.Sanity = T{:,:};

% get rmsd in tq, v and corr in delta;
Score = [];
for fq = fieldnames(FQ)'
	n=0;
	for my_norm = fieldnames(FQ.(fq{:}))'
		n=n+1;
		if strcmp(fq{:},'lfc')
			for g = 1:N_gene
				Score.(fq{:})(g,n) = corr(FQ.lfc.(my_norm{:})(g,:)',True.lfc(g,:)');
			end
		else
			Score.(fq{:})(:,n) = abs(FQ.(fq{:}).(my_norm{:}) - True.(fq{:}))./True.(fq{:}); 
		end
	end
end

% Get mean UMI per gene
T = readtable(['data/Simulated_Baron_Independent_Genes_UMI_count.txt'],'ReadRowNames',1,'delimiter','\t');
UMI = SC{:,:};
N_g = mean(UMI,2);

% Compute bins and find genes in bins
edges = [log10(0.0005) -2:1 log10(500)];
[n_per_bin,idx_bin] = histc(N_g,10.^edges);
n_per_bin = n_per_bin(1:end-1);
idx_bin = idx_bin(1:end-1);
N_bins = length(n_per_bin);

% Get score quantiled per bin
q = [.05 .25 .5 .75 .95];
Quant = [];
for b = 1:N_bins
	for fq = fieldnames(Score)'
		for n = 1:length(My_norm)
			for i = 1:length(q)
				Quant.(fq{:}).(My_norm{n})(i,b) = quantile(Score.(fq{:})(idx_bin==b,n),q(i));
			end
		end
	end
end

% plot paramters
tit.mean_tq = 'Errors in the estimates of mean expression';
ylab.mean_tq = {'relative error:' '|true-estimate|/true'};
tit.var = 'Errors in the estimates of expression variance';
ylab.var = {'relative error:' '|true-estimate|/true'};
tit.lfc = 'Correlations of true and estimated log fold-changes';
ylab.lfc = 'Pearson correlation';
xticklabel = {'0.0005' '0.01' '0.1' '1' '10' '500'};
my_ylim.mean_tq = [0.0004 36.7];
my_ylim.var = [0.0015 13.5];
my_ylim.lfc = [-0.04 1];

% sorting index
idx_sort.mean_tq = [3 10 9 8 5 1 2 7 4 6];
idx_sort.var =     [5 3 10 9 8 1 2 4 7 6];
idx_sort.lfc =     [5 3 10 8 9 1 2 4 7 6];

% bin position
bin_center = .5*(edges(1:end-1)+edges(2:end));
dx = linspace(-.4,.4,10);

% Plot box plots
figure('visible','off'); f=1;
for fq = fieldnames(Score)'

	subplot(2,2,f);f=f+1;
	hold on	
	k = 0;
	for n = idx_sort.(fq{:})
		k=k+1;
		x = bin_center + dx(k);

		% 5% - 95%
		plot([x;x],Quant.(fq{:}).(My_norm{n})([1 5],:),'linewidth',.5,'color',my_colors(n,:))
		% 25% - 75%
		plot([x;x],Quant.(fq{:}).(My_norm{n})([2 4],:),'linewidth',2.5,'color',my_colors(n,:))
		% 50%
		plot(x,Quant.(fq{:}).(My_norm{n})(3,:),'o','color',my_colors(n,:),'MarkerFaceColor',[1 1 1],'Markersize',2);
	end

	if ~strcmp(fq{:},'lfc')
		set(gca,'yscale','log')
		set(gca,'YTick',10.^[-3:1],'YTickLabel',{'0.001' '0.01' '0.1' '1' '10'})
	else
		set(gca,'YTick',0:.25:1)
	end
	set(gca,'XTick',edges,'XTickLabel',xticklabel);
	xlim([edges(1) edges(end)])
	ylim(my_ylim.(fq{:}))
	set(gca,'XGrid','on')
	xlabel('<UMI per cell>')
	ylabel(ylab.(fq{:}))
	title(tit.(fq{:}),'FontWeight','normal','FontSize',9)
	set(gca,'Fontsize',9);
end

% legend
y = [5:-1:1 5:-1:1];
x = [0*ones(1,5) 3*ones(1,5)];
subplot(2,2,4)
scatter(x,y,200,my_colors ,'s','filled')
hold on
for ii = 1:length(y)
	text(x(ii)+.5,y(ii),My_norm{ii},'FontSize',9)
end
axis([0 6 .5 6])
axis off

dim = [20 10];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['figure_6'],'-dpdf');
