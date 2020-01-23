clear all; close all; clc;

addpath('scripts')

My_norm = {'True','RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};

% define colors
my_colors = cbrewer('qual','Set1',length(My_norm)-1);
tmp_colors = cbrewer('qual','Dark2',8);
my_colors(6,:) = mean([my_colors(6,:); tmp_colors(6,:)]);
my_colors = [[0 0 0]; my_colors];

for n = 1:length(My_norm)
	if n==1
		% Need to run run_SimulatedBaron.m to create my_sim.mat
		load('data/SimulatedBaron/my_sim.mat');
		my_mean(:,n) = mean(E,2);
		my_cv(:,n) = std(E,0,2)./my_mean(:,n);
	else
		load(['data/SimulatedBaron_' My_norm{n} '_normalization_lin.mat']);
		M(isinf(M))=NaN;
		M(M<0) = 0;

		my_mean(:,n) = nanmean(M,2);
		my_cv(:,n) =  nanstd(M,0,2)./my_mean(:,n);
        % Normalize Sanity expression to the mean UMI count per cell
        if strcmp(My_norm{n},'Sanity')
            T = readtable(['data/SimulatedBaron_UMI_counts.txt'],'ReadRowNames',1,'delimiter', '\t');
            median_Nc = median(sum(T{:,:},1));
            M = median_Nc*M;
            my_mean(:,n) = my_mean(:,n)*median_Nc;
        end
	end
end

% CV-CV scatter
figure('visible','off');
for n = 2:length(My_norm)

	subplot(3,3,n-1);
	x = [.004 50];
	plot(x,x,'color',[.5 .5 .5]);
	hold on;
	scatter(my_cv(:,1),my_cv(:,n),1,log10(my_mean(:,1)),'.');
	set(gca,'xscale','log','yscale','log');
	grid on;	
	title([My_norm{n} ' ' num2str(corr(my_cv(:,1),my_cv(:,n)),2)],'FontWeight','Normal')
	axis([.004 50 .004 50])
	ticks = logspace(-2,1,4);
	set(gca,'XMinorGrid','off','YMinorGrid','off')
	set(gca,'ytick',ticks,'xtick',ticks);

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

dim = [20 17];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_2','-r1000','-dpng');

% Plot colorbar
figure;
sort_mean = sort(my_mean(:,1));
scatter(sort_mean([1:10 end]),sort_mean([1:10 end]),1,sort_mean([1:10 end]),'.')
hcb = colorbar;
set(get(hcb,'Title') ,'String','Gene mean expression');
hcb.Ruler.Scale = 'log';
hcb.Ruler.TickValues = [.001 .01 0.1 1 10 100 1000];

dim = [12 5];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_2_colorbar','-dsvg');

