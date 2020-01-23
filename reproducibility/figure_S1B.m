clear all; close all; clc;

addpath('scripts')

dataset = 'SimulatedBaron';
My_norm = {'True','RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};

% define colors
my_colors = cbrewer('qual','Set1',length(My_norm)-1);
tmp_colors = cbrewer('qual','Dark2',8);
my_colors(6,:) = mean([my_colors(6,:); tmp_colors(6,:)]);
my_colors = [[0 0 0]; my_colors];

my_mean = [];
my_cv = [];
for n = 1:length(My_norm)	
	if n==1
        % Need to run run_SimulatedBaron.m to create my_sim.mat
		load('data/SimulatedBaron/my_sim.mat');
		my_mean(:,n) = mean(E,2);
		my_cv(:,n) = std(E,0,2)./my_mean(:,n);
	else
		% Compute cv and mean
		load(['data/SimulatedBaron_' My_norm{n} '_normalization_lin.mat']);
		M(M<0) = 0;
		my_mean(:,n) = nanmean(M,2);
		my_cv(:,n) =  nanstd(M,0,2)./my_mean(n,:);

		% Normalize Sanity mean expression to the mean UMI count per cell
		if strcmp(My_norm{n},'Sanity')
			T = readtable('data/SimulatedBaron_UMI_counts.txt','ReadRowNames',1,'delimiter', '\t');
			median_Nc = median(sum(T{:,:},1));
			my_mean(:,n) = my_mean(:,n)*median_Nc;
		end
	end
end

% Plot scatter cv-mean
figure('visible','off');
for n = 1:length(My_norm)
	fprintf([My_norm{n} '\n'])
	subplot(1,length(My_norm),n);
	hold on;
	scatter(my_mean(:,n),my_cv(:,n),1,my_colors(n,:),'.');
	set(gca,'xscale','log','yscale','log');
	grid off;
	title({My_norm{n},num2str(corr(log2(my_mean(:,n)),log2(my_cv(:,n))),2)},'FontWeight','Normal')

	if strcmp(My_norm{n},'scVI')
		axis([1e7 4e14 .004 50])
		xtick = logspace(8,14,3);
	else
		axis([1e-4 1e4 0.004 50])
		xtick = logspace(-3,3,3);
	end
	set(gca,'xtick',xtick)

	if n==1
		ylabel({'Simulated','CV'})
		ytick = logspace(-2,1,4);
		set(gca,'ytick',ytick)
	else
		set(gca,'ytick',[])
	end

	if n == 5
		xlabel('mean expression level')
	end
end

dim = [43 5];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S1B'],'-dsvg');
