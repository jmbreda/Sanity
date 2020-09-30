clear all; close all;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','scVI'};
my_colors = load('data/my_colors.txt');
my_colors(9,:) = [];

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo' 'Simulated_Baron_Independent_Genes'};
Datasets_label = {'Zeisel' 'Baron' 'Chen' 'LaManno Embryo' 'LaManno ES' 'LaManno MouseEmbryo' 'Simulated'};

figure('visible','off');
k=0;
for d = 1:length(Datasets)
	for n = 1:length(My_norm)
		k=k+1;
		
		% Loag gene expression
		load(['data/' Datasets{d} '_' My_norm{n} '_normalization_lin.mat']);

		% Compute CV and mean
		M(M<0) = 0;
		my_mean = nanmean(M,2);
		my_cv =  nanstd(M,0,2)./my_mean;

		% plot cv mean scatter
		subplot(length(Datasets),length(My_norm),k);
		scatter(my_mean,my_cv,1,my_colors(n,:),'.');
		set(gca,'xscale','log','yscale','log');
		grid off;
		
		if d==1
			title({My_norm{n},num2str(real(corr(log(my_mean),log(my_cv))),2)},'FontWeight','Normal')
		else
			title(num2str(real(corr(log(my_mean),log(my_cv))),2),'FontWeight','Normal')
		end
		ylim([0.005 125])
		if n==1
			ylabel({Datasets_label{d},'CV'})
			ytick = logspace(-2,2,3);
			yticklabel = {'0.01', '1','100'};
			set(gca,'ytick',ytick,'yticklabel',yticklabel);
		else
			set(gca,'ytick',ytick,'yticklabel',[])
		end

		if strcmp(My_norm{n},'scVI') || strcmp(My_norm{n},'Sanity') || strcmp(My_norm{n},'True')
			xlim([1e-8 1e-2])
			xtick = logspace(-8,-2,3);
		else	
			xlim([1e-3 1e2])
			xtick = logspace(-2,2,3);
		end

		if d == length(Datasets)
			set(gca,'xtick',xtick)
			if n==5
				xlabel('mean expression level')
			end
		else
			set(gca,'xticklabel',[])
		end
	end
end
dim = [40 30];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S2'],'-dpng');

% Print for simualted data only
% Need to run run_Simulations.m to create data/Simulated_Baron_Independent_Genes.mat
load('data/Simulated_Baron_Independent_Genes.mat');
my_mean = mean(exp(E),2);
my_cv = std(exp(E),0,2)./my_mean;

% plot cv mean scatter
figure('visible','off');
scatter(my_mean,my_cv,1,'k.');
set(gca,'xscale','log','yscale','log');
grid off;
title({'True',num2str(corr(log(my_mean),log(my_cv)),2)},'FontWeight','Normal')
ylim([0.005 125])
ylabel({'Simulated','CV'})
ytick = logspace(-2,2,3);
yticklabel = {'0.01', '1','100'};
set(gca,'ytick',ytick,'yticklabel',yticklabel);
xlim([1e-8 1e-2])
xtick = logspace(-8,-2,3);
set(gca,'xtick',xtick)
xlabel('mean expression level')
dim = [6 6];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S2_True'],'-dpng');
