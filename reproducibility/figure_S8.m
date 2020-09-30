clear all; close all;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};

% Get true correlation matrix in simulation 
load('data/Simulated_Baron_Independent_Genes.mat');
C = corr(E');
rho{1} = C;
clear c;
c{1} = rho{1} - 3*tril(ones(size(rho{1})));
c{1} = c{1}(c{1}>=-1);

% Compare with inferred correlation matrix
figure('visible','off'); f=1;
for meth_2 = My_norm
	
	% Load gene expression and compute et correlation matrix
	load(['data/' Datasets{k} '_' My_norm{n} '_normalization.mat']);
	C = corr(M');
	rho{2} = C;
	c{2} = rho{2} - 3*tril(ones(size(rho{2})));
	c{2} = c{2}(c{2}>=-1);

	% Create 2d histogram
	Edges = [-1.005:0.01:1.005]';
	[H X] = hist3([c{2} c{1}],'Edges',{Edges Edges});

	% Plot
	subplot(4,3,f);
	f=f+1;
	imagesc(X{2},X{1},log10(H),'AlphaData',~(H==0))
	my_map = cbrewer('seq','Blues',256);
	colormap(my_map(64:end,:))
	axis xy;
	hcb = colorbar;
	my_ticks = 0:2:floor( max(max(log10(H))) );
	for i = 1:length(my_ticks)
		my_ticklabels{i} = ['10^' num2str(my_ticks(i))];
	end
	set(hcb,'Ticks',my_ticks,'TickLabels',my_ticklabels);
	set(get(hcb,'Title') ,'String',{'nr. of' 'gene pairs'});
	hold on;
	plot([-1 0; 1 0],[0 -1; 0 1],'r:')
	axis([-1 1 -1 1])
	title(meth_2{:},'FontWeight','normal')
    if f==9
        xlabel(['True Pearson correlation'])
    end
    if f==5
        ylabel(['Inferred Pearson correlation'])
    end
end

dim = [20 20];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S8','-dpdf');

