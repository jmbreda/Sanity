clear all; close all;

addpath('scripts')

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo' 'SimulatedBaron'};
My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};
Subplot_cols = length(My_norm)+1;

% define colors
tmp_colors = cbrewer('qual','Set1',length(My_norm));
for ii = 1:length(My_norm)
    my_colors{ii} = tmp_colors(ii,:);
end
tmp_colors = cbrewer('qual','Dark2',8);
my_colors{6} = mean([my_colors{6}; tmp_colors(6,:)]);
my_colors{9} = tmp_colors(8,:);

figure('visible','off')

for d = 1:length(Datasets)

	% Get total UMI count per cell
	T = readtable(['data/' Datasets{d} '_UMI_counts.txt'],'ReadRowNames',1,'delimiter', '\t');
	UMI_per_cell = sum(T{:,:},1)';

	% Mesure correlation on various normalisations
	C = nan(size(UMI_per_cell,1),length(My_norm));
	for n = 1:length(My_norm)
		load(['data/' Datasets{d} '_' My_norm{n} '_normalization.mat'])
		C(:,n) = corr(M',log(UMI_per_cell));
	end

	x_names = [];

	if strcmp(Datasets{d},'SimulatedBaron') 
		idx_subplot = Subplot_cols*(d-1)+1:d*Subplot_cols;

		% Need to run run_SimulatedBaron.m to create my_sim.mat
		load('data/SimulatedBaron/my_sim.mat');
		C = [corr(E',N_c') C];
	else
		idx_subplot = Subplot_cols*(d-1)+2:d*Subplot_cols;
	end

	subplot(length(Datasets),Subplot_cols,idx_subplot)
	
	if strcmp(Datasets{d},'SimulatedBaron')
		distributionPlot(C,'color',[[0 0 0] my_colors],'showMM',0,'globalNorm',1,'distWidth',1);
		axis([0.5 length(My_norm)+1.5 -1 1])
	else
		distributionPlot(C,'color',my_colors,'showMM',0,'globalNorm',1,'distWidth',1);
		axis([0.5 length(My_norm)+.5 -1 1])
	end
	grid on

	if d < length(Datasets)
		set(gca,'Xtick',1:length(My_norm),'XTickLabel',[])
	end
	
	if d==4
		set(get(gca,'Ylabel'),'String','Correlation between gene log expression and log(total UMI count per cell)')
	end
	title(Datasets{d},'FontWeight','normal')
end

set(gca,'XTickLabel',['True' My_norm],'XTickLabelRotation',0)

dim = [24 28];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_S2','-dpdf');

