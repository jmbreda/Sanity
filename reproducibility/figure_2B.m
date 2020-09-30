clear all; close all;

addpath('scripts')
My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','scVI'};
Datasets = {'SC_2i','SC_serum','RNA_2i','RNA_serum'};

% define colors
C = lines(7);
C = C([1 6 7 2],:);

% Compute mean and CV
for d = 1:length(Datasets)
	for n = 1:length(My_norm)
		load(['data/Gruen_ESC_' Datasets{d} '_' My_norm{n} '_normalization_lin.mat']);

        if strcmp(My_norm{n},'sctransform')
			my_mean = readtable(['data/Gruen_SC_2i_sctransform_mean_expression.txt']);
			MU.(My_norm{n}).(Datasets{d}) = my_mean{:};
		else
			MU.(My_norm{n}).(Datasets{d}) = nanmean(M,2);
		end
		CV.(My_norm{n}).(Datasets{d}) = nanstd(M,0,2)./MU.(My_norm{n}).(Datasets{d});
	end
end

% Scatter plot
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
	title([My_norm{n} ' Corr: ' num2str(rho,2)],'FontWeight','Normal')

	if strcmp(My_norm{n},'Sanity')
		xlim([4e-7 1e-2])
		xtick = logspace(-6,-2,3);
		xticklabel = {'0.000001' '0.0001' '0.01'};
	elseif strcmp(My_norm{n},'scVI')
		xlim([1e-6 1e-2])
		xtick = logspace(-6,-2,3);
		xticklabel = {'0.000001' '0.0001' '0.01'};
	elseif strcmp(My_norm{n},'sctransform')
		xlim([1e-2 4e5])
		xtick = logspace(-1,5,3);
		xticklabel = {'0.1' '100' '10000'};
	else
		xlim([1e-2 1e3])
		xtick = logspace(-1,3,3);
		xticklabel = {'0.1' '10' '1000'};
	end
	set(gca,'xtick',xtick,'xticklabel',xticklabel)

	ytick = logspace(-2,1,4);
	yticklabel = {'0.01' '0.1' '1' '10'};
	set(gca,'ytick',ytick,'yticklabel',yticklabel);
	if  mod(n,3)~=1
		set(gca,'yticklabel',[])
	end

	if n==8
		xlabel('mean expression')
	end
	if n==4
		ylabel('CV')
	end
end
dim = [24 18];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_2B','-dpdf');
