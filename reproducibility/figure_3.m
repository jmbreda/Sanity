clear all; close all; clc;

addpath('scripts')

My_norm = {'True' 'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute',            'sctransform','scVI'};

for n = 1:length(My_norm)
	if strcmp(My_norm{n},'True')
		% Need to run run_Simulations.m to create data/Simulated_Baron_Independent_Genes.mat
		load('data/Simulated_Baron_Independent_Genes.mat');
		my_var(:,n) =  var(E,0,2);
	elseif strcmp(My_norm{n},'Sanity')
		tmp = readtable('data/Simulated_Baron_Independent_Genes_Sanity_variance.txt');
		my_var(:,n) = tmp{:,:};
	else
		load(['data/Simulated_Baron_Independent_Genes_' My_norm{n} '_normalization.mat']);
		my_var(:,n) =  nanvar(M,0,2);
	end
end

% var-var scatter
figure('visible','off');
for n = 2:length(My_norm)
	subplot(4,3,n-1);
	x = [0 24];
	plot(x,x,'color',[.5 .5 .5]);
	hold on;
	scatter(my_var(:,1),my_var(:,n),1,log10(mean_UMI),'.');
	if strcmp(My_norm{n},'Deconvolution')
		title(['Deconv. r: ' num2str(corr(my_var(:,1),my_var(:,n)),1)],'FontWeight','Normal')
	else
		title([My_norm{n} ' r: ' num2str(corr(my_var(:,1),my_var(:,n)),2)],'FontWeight','Normal')
	end

	x_min = 0;
	x_max = 24;
	y_min = 0;
	y_max = 11;
	axis([x_min x_max y_min y_max])
	xtick = 0:5:20;
	ytick = 0:2.5:10;
	set(gca,'ytick',ytick,'xtick',xtick);

	if strcmp(My_norm{n},'sctransform')
		axis([x_min x_max y_min max(my_var(:,n))]);
		set(gca,'ytick',0:15:60);
	elseif mod(n-1,3)~=1
		set(gca,'yticklabel',[]);
	end

	if n-1<7
		set(gca,'xticklabel',[]);
	end

	if n-1==4
		ylabel('Inferred variance')
	end
	if n-1==8
		xlabel('True variance')
	end
end
dim = [16 20];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_3','-r1000','-dpng');

% Plot colorbar
figure;
sort_mean = sort(my_mean(:,1));
scatter(sort_mean([1:10 end]),sort_mean([1:10 end]),1,sort_mean([1:10 end]),'.')
hcb = colorbar;
set(get(hcb,'Title') ,'String','<UMI> per cell');
hcb.Ruler.Scale = 'log';
hcb.Ruler.TickValues = [.001 .01 0.1 1 10 100 1000];

dim = [12 5];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_3_colorbar','-dsvg');

