clear all; close all; clc;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','SanityErrorbar','SAVER','scImpute','sctransform','scVI' 'Sanity'};

% Get true distance from simulated data (run run_Simulations.m)
load('data/Simulated_Branched_Random_Walk.mat');
Eucl{1} = D_true;
Edges_1 = linspace(0,max(Eucl{1}),200)';

% Compare with inferred distance
figure('visible','off'); f=1;
for my_norm = My_norm

	% Compute cell-cell distances
	if strcmp(my_norm{:},'SanityErrorbar')
		% For Sanity distance with errorbar, run:
		% Sanity_distance -f 'Sanity output folder' -s2n 0 -err 1
		Eucl{2} = load('path/to/Sanity/folder/cell_cell_distance_with_errorbar.txt');
	else
		load(['data/Simulated_Branched_Random_Walk_' my_norm{:} '_normalization.mat']);
		Eucl{2} = pdist(M');
	end

	% Compute 2d histogram
	Edges_2 = linspace(0,max(Eucl{2}),200)';
	[H X] = hist3([Eucl{2} Eucl{1}],'Edges',{Edges_2 Edges_1});

	% Plot
	subplot(4,3,f);
	f=f+1;
	x = [0 max(Edges_2)];
	plot(x,x,'color',[.5 .5 .5])
	hold on
	imagesc(X{2},X{1},log10(H),'AlphaData',~(H==0))
	my_map = cbrewer('seq','Blues',256);
	colormap(my_map(64:end,:))
	axis xy;
	hcb = colorbar;
	my_ticks = 0:3;
	my_ticklabels = {'1' '10' '100' '1000'};
	set(hcb,'Ticks',my_ticks,'TickLabels',my_ticklabels);
	set(get(hcb,'Title') ,'String',{'nr. of' 'cell pairs'});
	axis([0 Edges_1(end) 0 Edges_2(end)]);
	xlabel('True Euclidean distance')
	ylabel('Inferred Euclidean distance')

	if strcmp(my_norm{:},'SanityErrorbar')
		name = 'Sanity (with error bars) Corr: ';
	elseif strcmp(my_norm{:},'Sanity')
		name = 'Sanity (without error bars) Corr: ';
	else
		name = [my_norm{:} ' Corr: '];
	end
	title([name num2str(corr(Eucl{1},Eucl{2}),2)],'FontWeight','normal')
end

% Print
dim = [36 36];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S10','-dpdf');

