clear all; close all; clc;

addpath('scripts')
addpath('scripts/bhtsne')

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo'};
fig_name = {'5A' 'S13' 'S14' 'S15' 'S16' 'S17'};
My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};

% Parameters for tsne
numDims = 2;
pcaDims = 50;
theta = .5;
alg = 'svd';

for k = 1:length(Datasets)
	% Get celltype annotation
	Celltype = textread(['data/' Datasets{k} '_Celltype.txt'],'%s\n');
	[Celltype,~,G] = unique(Celltype);

	% define colors
	if k==1
		my_colors = distinguishable_colors(max(G)+7);
		my_colors = my_colors(8:end,:);
	else
		my_colors = distinguishable_colors(max(G));
	end

	figure('visible','off'); f=1;
	% Get tsne transformation and plot
	for my_norm = My_norm;
		
		% Compute tsne on normalization
		load(['data/' Datasets{k} '_' my_norm{:} '_normalization.mat']);
		N_cells = size(M,2);
		N_classes = length(Celltype);
		perplexity = N_cells/N_classes;	
		tsne_map = fast_tsne(M', numDims, pcaDims, perplexity, theta, alg);
		
		% Plot tsne representation
		subplot(2,6,f);f=f+1;
		if f==6
			f=f+1;
		end
		gscatter(tsne_map(:,1),tsne_map(:,2),G,my_colors,'.','3')
		legend off;
		box off
		set(gca,'xtick',[],'ytick',[])
		title(my_norm{:},'interpreter','none','FontWeight','normal')
		axis([min(tsne_map(:,1)) max(tsne_map(:,1)) min(tsne_map(:,2)) max(tsne_map(:,2))])
	end
	
	% Plot legend
	subplot(2,6,[6 12])
	y = length(Celltype):-1:1;
	x = zeros(size(y));
	x_lim = [-.5 5];
	switch k
	case 1
		y_lim = [-3 11];
	case 2
		y_lim = [-1 16];
	case 3
		y_lim = [-1 16];
	otherwise
		y_lim = [min(y)-1 max(y)+1];
	end
	scatter(x,y,60,my_colors,'s','filled')
	hold on
	for ii = 1:length(y)
		text(x(ii)+.5,y(ii),strrep(Celltype{ii},'_',' '),'FontSize',8)
	end
	xlim(x_lim)
	ylim(y_lim)
	axis off
	
	suptitle(Datasets{k})

	dim=[29 10];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim])
	print(gcf,['figure_' fig_name{k} ],'-dpdf')
end
