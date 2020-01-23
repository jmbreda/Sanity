clear all; close all;

addpath('scripts');

My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};
Datasets = {'SC_2i','SC_serum','RNA_2i','RNA_serum'};
Datasets_label = {'Single-cell 2i','Single-cell serum','Aliquots 2i','Aliquots serum'};

% define colors
C = lines(7);
C = C([1 6 7 2],:);
dark_blue = C(1,:);
light_blue = C(2,:);
dark_red = C(3,:);
light_red = C(4,:);
my_color.SC_2i = dark_blue;
my_color.SC_serum = light_blue;
my_color.RNA_2i = dark_red;
my_color.RNA_serum = light_red;

% Compute CV and quantiles at .05 .25 .5 .75 .95
for d = 1:length(Datasets)
	for n = 1:length(My_norm)
		load(['data/Gruen_ESC_' Datasets{d} '_' My_norm{n} '_normalization_lin.mat']);

		% Normalize Sanity expression to the mean UMI count per cell
		if strcmp(My_norm{n},'Sanity')
			T = readtable(['data/Gruen_ESC_' Datasets{d} '_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
			mean_Nc = mean(sum(T{:,:},1));
			M = mean_Nc*M;
		end
		
		CV.(Datasets{d})(:,n) = nanstd(M,0,2)./nanmean(M,2);
	end
	
	Quant.(Datasets{d}) = zeros(5,length(My_norm));
	q = [.05 .25 .5 .75 .95];
	for i = 1:5
		Quant.(Datasets{d})(i,:) = quantile(CV.(Datasets{d}),q(i)); 
	end
end

% Boxplot
xticks = 1:length(My_norm);
dx = linspace(-0.25,0.25,4);
clear x;
x.SC_2i = xticks + dx(1);
x.SC_serum = xticks + dx(2);
x.RNA_2i = xticks + dx(3);
x.RNA_serum = xticks + dx(4);

figure('visible','off')
hold on
for d = Datasets
	% 5% - 95%
	plot([x.(d{:});x.(d{:})],Quant.(d{:})([1 5],:),'linewidth',1,'color',my_color.(d{:}))
	% 25% - 75%
	plot([x.(d{:});x.(d{:})],Quant.(d{:})([2 4],:),'linewidth',8,'color',my_color.(d{:}))
	% 50%
	plot(x.(d{:}),Quant.(d{:})(3,:),'o','color',my_color.(d{:}),'MarkerFaceColor',[1 1 1])
	plot(x.(d{:}),Quant.(d{:})(3,:),'.','color',my_color.(d{:}))
end
set(gca,'yscale','log','ytick',[.05 .1 .2 .5 1 2 5 10],'yticklabel',{'0.05' '0.1' '0.2' '0.5' '1' '2'  '5' '10'})
axis([0.5 length(My_norm)+.5 0.03 10])
set(gca,'XTick',xticks,'XTickLabel',My_norm)
set(gca,'Ygrid','on')
box on
ylabel('CV')

dim = [24 10];  
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_1A','-dpdf');

