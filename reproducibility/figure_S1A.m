clear all; close all;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};
Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo'};
Datasets_label = {'Zeisel' 'Baron' 'Chen' 'LaManno Embryo' 'LaManno ES' 'LaManno MouseEmbryo' 'Simulated'};


% define colors
my_colors = cbrewer('qual','Set1',length(My_norm));
tmp_colors = cbrewer('qual','Dark2',8);
my_colors(6,:) = mean([my_colors(6,:); tmp_colors(6,:)]);

figure('visible','off');
k=0;
for d = 1:length(Datasets)
	for n = 1:length(My_norm) 
		k=k+1;
		
		% Compute cv and mean
		load(['data/' Datasets{d} '_' My_norm{n} '_normalization_lin.mat']);
		M(M<0) = 0;
		my_mean = nanmean(M,2);
		my_cv =  nanstd(M,0,2)./my_mean;

	    % Normalize Sanity mean expression to the mean UMI count per cell
        if strcmp(My_norm{n},'Sanity')
            T = readtable(['data/' Datasets{d} '_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
			median_Nc = median(sum(T{:,:},1));
			my_mean = my_mean*median_Nc;
		end

		% plot cv mean scatter
		subplot(length(Datasets),length(My_norm),k);
		hold on;
		scatter(my_mean,my_cv,1,my_colors(n,:),'.');
		set(gca,'xscale','log','yscale','log');
		grid off;
		
		if d==1
			title({My_norm{n},num2str(corr(log2(my_mean),log2(my_cv)),2)},'FontWeight','Normal')
		else
			title(num2str(corr(log2(my_mean),log2(my_cv)),2),'FontWeight','Normal')
		end

		if ~strcmp(My_norm{n},'scVI')
			axis([1e-3 1e2 0.005 120.2])
			xtick = logspace(-2,2,3);

			if n==1
				ylabel({Datasets_label{d},'CV'})

				ytick = logspace(-2,2,5);
				set(gca,'ytick',ytick)
			else
				set(gca,'ytick',[])
			end

			if d == length(Datasets)
				set(gca,'xtick',xtick)
				if n==4
					xlabel('mean expression level')
				end
			else
				set(gca,'xtick',[])
			end
		else
			axis([min(my_mean) max(my_mean) min(my_cv) max(my_cv)])
			axis([min(my_mean) max(my_mean) 0.005 2e2])

			xtick = logspace(-2,10,7);
			if d==3
				xtick = logspace(-2,6,3);
			elseif d==5
				xtick = [1e4 1e7 1e10];
			elseif d==7
				xtick = [1e8 1e11 1e14];
			end
			set(gca,'xtick',xtick)

			set(gca,'ytick',[])

			if 0 && d==7 && n==1
				axis([1e-3 1e4 0.005 2e2])
			end
		end
	end
end

dim = [40 30];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S1A'],'-dpdf');


