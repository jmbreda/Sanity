clear all; close all; clc;

addpath('scripts')

My_norm = {'True','RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};
dataset = 'Simulated_Baron_Independent_Genes';

for n = 1:length(My_norm)
	if n==1
		% Need to run run_Simulations.m to create data/Simulated_Baron_Independent_Genes.mat
		load(['data/' dataset '.mat']);
		my_mean(:,n) = mean(exp(E),2);
		my_cv(:,n) = std(exp(E),0,2)./my_mean(:,n);
	
		% Compute mean UMI count per cell
		mean_UMI = mean(UMI,2);
	else
		% Load gene expression and compute mean and cv
		load(['data/' dataset '_' My_norm{n} '_normalization_lin.mat']);
		M(isinf(M))=NaN;
		M(M<0) = 0;
		my_mean(:,n) = nanmean(M,2);
		my_cv(:,n) =  nanstd(M,0,2)./my_mean(:,n);
	end
end

% CV-CV scatter
figure('visible','off');
for n = 2:length(My_norm)

	subplot(3,3,n-1);
	x = [.004 230];
	plot(x,x,'color',[.5 .5 .5]);
	hold on;
	scatter(my_cv(:,1),my_cv(:,n),1,log10(mean_UMI),'.');
	set(gca,'xscale','log','yscale','log');
	grid on;	
	title([My_norm{n} ' Corr: ' num2str(corr(my_cv(:,1),my_cv(:,n)),2)],'FontWeight','Normal')
	axis([.002 230 .002 230])
	ticks = logspace(-2,2,5);
	ticklabels = {'0.01' '0.1' '1' '10' '100'};
	set(gca,'XMinorGrid','off','YMinorGrid','off')
	set(gca,'ytick',ticks,'yticklabel',ticklabels,'xtick',ticks,'xticklabel',ticklabels);

	if mod(n-1,3)~=1
		set(gca,'yticklabel',[]);
	end
	if n-1<7
		set(gca,'xticklabel',[]);
	end
	if n-1==4
		ylabel('Inferred CV')
	end
	if n-1==8
		xlabel('True CV')
	end
end
dim = [16 15];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S3','-r1000','-dpng');

% Plot colorbar
figure;
sort_mean = sort(mean_UMI);
scatter(sort_mean([1:10 end]),sort_mean([1:10 end]),1,sort_mean([1:10 end]),'.')
hcb = colorbar;
set(get(hcb,'Title') ,'String','<UMI> per cell');
hcb.Ruler.Scale = 'log';
hcb.Ruler.TickValues = [.001 .01 0.1 1 10 100 1000];

dim = [12 5];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_S3_colorbar','-dsvg');

