clear all; close all;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};

% define colors
tmp_colors = load('data/my_colors.txt');
for ii = 1:length(My_norm)
	my_colors{ii} = tmp_colors(ii,:);
end

% Get total UMI count per cell
T = readtable(['data/Zeisel_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
UMI_per_cell = sum(T{:,:},1)';

figure('visible','off')

% Example gene : Zbed3 in Zeisel
idx_Zbed3 = 8709;
for n = 1:length(My_norm)

	load(['data/Zeisel_' My_norm{n} '_normalization.mat'])
	rho = corr(M(idx_Zbed3,:)',log(UMI_per_cell));

	subplot(4,length(My_norm),n)
	semilogx(UMI_per_cell,M(idx_Zbed3,:),'.','color',my_colors{n},'MarkerSize',3)
	title(['R = ' num2str(rho,2)],'FontWeight','normal','Fontsize',8)
	set(gca,'xtick',[],'ytick',[])
	box off
	axis([min(UMI_per_cell) max(UMI_per_cell) min(M(idx_Zbed3,:)) max(M(idx_Zbed3,:))])
	if n==5
		xlabel('log N_c')
	end
	if n==1
		ylabel({'Zbed3','expression'})
	end
end

% Get correlation between gene expression and total UMI counts

% Mesure correlation on the different normalisations
C = nan(size(UMI_per_cell,1),length(My_norm));
for n = 1:length(My_norm)
	load(['data/Zeisel_' My_norm{n} '_normalization.mat'])
	C(:,n) = corr(M',log(UMI_per_cell));
end

% plot violin distributions
subplot(4,1,[2 4])
distributionPlot(C,'color',my_colors,'showMM',0,'globalNorm',1,'distWidth',1);
grid on
set(gca,'Xtick',1:length(My_norm),'XTickLabel',[])
axis([0.5 length(My_norm)+.5 -1 1])
set(get(gca,'Ylabel'),'String',{'Correlation between','gene log expression','and log total UMI count'})
set(gca,'XTickLabel',My_norm)
dim = [24 10];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_4'],'-dpdf');
