clear all; close all;

addpath('scripts')

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo' 'SimulatedBaron'};
My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};

% define colors
my_colors = cbrewer('qual','Set1',length(My_norm));
tmp_colors = cbrewer('qual','Dark2',8);
my_colors(6,:) = mean([my_colors(6,:); tmp_colors(6,:)]);
my_colors = my_colors;
set(groot,'defaultAxesColorOrder',my_colors)

% Plot distribution of pairwise correlation
figure('visible','off')
for k = 1:length(Datasets)
	
	H = [];
	X = [];
	for n = 1:length(My_norm);
		% Compute gene pairwise correlation
		load(['data/' Datasets{k} '_' My_norm{n} '_normalization.mat']);
		C = corr(M');

		% Get distribution
		idx = find( triu(ones(size(C)),1) );
		[h,x] = hist(C(idx),200);
		save(out_file,'h','x');
		H(:,n) = h';
		X(:,n) = x';
	end

	% Normalize distribution
	D = bsxfun(@rdivide,H,sum(H,1));

	% Plot distribution 
	subplot(3,3,k);
	semilogy(X,D,'linewidth',1.5);
	if k>3
		xlabel('\rho')
	else
		set(gca,'XtickLabel',[]);
	end
	if mod(k-1,3)==0
		ylabel('density')
		set(gca,'Ytick',[1e-8 1e-6 1e-4 1e-2 1e0])
	else
		set(gca,'YtickLabel',[]);
	end
	title(Datasets{k},'FontWeight','normal')
	axis([-1 1 1e-6 1])
end

% plot Legend :
y = length(My_norm):-1:1;
x = 1*ones(size(y));

subplot(3,3,9)
scatter(x,y,60,my_colors,'s','filled')
xlim([0 3])
hold on
for ii = 1:length(y)
	text(x(ii)+.2,y(ii),My_norm{ii},'FontSize',12)
end
axis off

dim = [32 24];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S4'],'-dpdf');
