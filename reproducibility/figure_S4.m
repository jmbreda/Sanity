clear all; close all;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform',     'scVI'};
Subplot_cols = length(My_norm)+1;

% define colors
tmp_colors = load('data/my_colors.txt');
for ii = 1:length(My_norm)
	my_colors{ii} = tmp_colors(ii,:);
end

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo' 'Simulated_Baron_Independent_Genes'};
Datasets_name = {'Zeisel' 'Baron' 'Chen' 'LaManno Embryo' 'LaManno ES' 'LaManno MouseEmbryo' 'Simulated'};

figure('visible','off')
for d = 1:length(Datasets)

	% Add "true" if simulated dataset
	sim = strcmp(Datasets{d},'Simulated_Baron_Independent_Genes');
	if sim
		My_norm = ['True' My_norm];
	end

	% Get total UMI count per cell
	T = readtable(['data/' Datasets{d} '_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
	UMI_per_cell = sum(T{:,:},1)';

    % Mesure correlation on various normalisations
    C = nan(size(UMI_per_cell,1),length(My_norm));
    for n = 1:length(My_norm)
		if sim & n==1
			load(['data/Simulated_Baron_Independent_Genes.mat']);
			C(:,n) = corr(E',log(UMI_per_cell));
		else
			load(['data/' Datasets{d} '_' My_norm{n} '_normalization.mat'])
			C(:,n) = corr(M',log(UMI_per_cell));
		end
    end

	% Get subplot indices
	if sim
		idx_subplot = Subplot_cols*(d-1)+1:d*Subplot_cols;
	else
		idx_subplot = Subplot_cols*(d-1)+2:d*Subplot_cols;
	end
	subplot(length(Datasets),Subplot_cols,idx_subplot)
	
	% plot violin distributions
	if sim
		distributionPlot(C,'color',[[0 0 0] my_colors],'showMM',0,'globalNorm',1,'distWidth',1);
	else
		distributionPlot(C,'color',my_colors,'showMM',0,'globalNorm',1,'distWidth',1);
	end
	axis([0.5 length(My_norm)+.5 -1 1])
	grid on
	if d < length(Datasets)
		set(gca,'Xtick',1:length(My_norm),'XTickLabel',[])
	end
	if d==4
		set(get(gca,'Ylabel'),'String','Correlation between gene log expression and log(total UMI count per cell)')
	end
	title(Datasets_name{d},'FontWeight','normal')
end
set(gca,'XTickLabel',My_norm)
dim = [24 28];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S4'],'-dpdf');
