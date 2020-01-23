clear all; close all;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};
Datasets = {'SC_2i','SC_serum','RNA_2i','RNA_serum'};

% define colors
C = lines(7);
C = C([1 6 7 2],:);

% Compute mean and CV
for d = 1:length(Datasets)
	for n = 1:length(My_norm)
		fprintf([Datasets{d} ' ' My_norm{n} '\n'])

		load(['data/Gruen_ESC_' Datasets{d} '_' My_norm{n} '_normalization_lin.mat']);

		% Normalize Sanity expression to the mean UMI count per cell
		if strcmp(My_norm{n},'Sanity')
			T = readtable(['data/Gruen_ESC_' Datasets{d} '_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
			mean_Nc = mean(sum(T{:,:},1));
			M = mean_Nc*M;
		end

		MU.(My_norm{n}).(Datasets{d}) = nanmean(M,2);
		CV.(My_norm{n}).(Datasets{d}) = nanstd(M,0,2)./MU.(My_norm{n}).(Datasets{d});
	end
end

% Scatter 3x3
figure('visible','off');
for n = 1:length(My_norm)
	my_mu = [];
	my_cv = [];
	G = [];
	for d = 1:length(Datasets)
		my_mu = [my_mu; MU.(My_norm{n}).(Datasets{d})];
		my_cv = [my_cv; CV.(My_norm{n}).(Datasets{d})];
		G = [G; d*ones(length(MU.(My_norm{n}).(Datasets{d})),1)];
	end

	% Random permutation for colors
	idx_perm = randperm(length(G));

	subplot(3,3,n);
	scatter(my_mu(idx_perm),my_cv(idx_perm),1,C(G(idx_perm),:),'.');
	set(gca,'xscale','log','yscale','log');
	grid off;

	rho = corr(log(my_mu),log(my_cv));
	title([My_norm{n} ' ' num2str(rho,2)],'FontWeight','Normal')

	if strcmp(My_norm{n},'scVI')
		axis([1e0 1e18 0.005 10])
		xtick = [1e0 1e9 1e18];
	else
		axis([0.01 1e3 0.005 10])
		xtick = logspace(-1,3,3);
		xtick = [1e-1 1e1 1e3];
	end
	set(gca,'xtick',xtick)
	ytick = logspace(-2,1,4);
	set(gca,'ytick',ytick)

	if  mod(n,3)==1
		set(gca,'ytick',ytick)
	else
		set(gca,'ytick',[])
	end

	if n==8
		xlabel('mean expression level')
	end
	if n==4
		ylabel('CV')
	end
end
dim = [24 18];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_1B','-dpdf');

