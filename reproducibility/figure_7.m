clear all; close all; clc;

addpath('~/Matlab_Scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI' 'SanityErrorbar'};
load('data/my_colors.txt');
K = [2:2:10 15:5:50];

% get knn
load('data/Simulated_Branched_Random_Walk.mat')% from run_Simulations.m
[~,my_knn] = sort(squareform(D));
knn.True = my_knn(2:end);
for my_norm = My_norm
	if strcmp(My_norm{:},'SanityErrorbar')
		% For Sanity distance with errorbar, run:
		% Sanity_distance -f 'Sanity output folder' -s2n 0 -err 1
		my_D = load('path/to/Sanity/folder/cell_cell_distance_with_errorbar.txt');
	else
		load(['data/Simulated_Branched_Random_Walk_' my_norm{:} '_normalization.mat']);
		my_D = pdist(M');
	end
	[~,my_knn] = sort(squareform(my_D));
	knn.(my_norm{:}) = my_knn(2:end,:);
end

% Compute fraction of correct nearest neighbors
N_cell = size(knn.True,2);
f_nn_mean = zeros(length(K),length(My_norm));
for jj = 1:(length(My_norm));
	for ii = 1:length(K)
		f_nn_mean(ii,jj) = 0;
		for c = 1:N_cell
			N_nn = length(intersect(knn.True(1:K(ii),c),knn.(My_norm{jj})(1:K(ii),c)));
			f_nn_mean(ii,jj) = f_nn_mean(ii,jj) + N_nn/K(ii);
		end
		f_nn_mean(ii,jj) = f_nn_mean(ii,jj)/N_cell;
	end
end

% Plot fraction of correct nearest neighbors
figure('visible','off')
subplot(1,3,1)
plot(K',f_nn_mean,'.-');
grid on
set(gca,'ColorOrder',my_colors);
axis([0 max(K) 0 1])
xlabel('Number of nearest neighbors')
ylabel('Fraction of correct neighbors')
title('All genes','FontWeight','normal')

% Compute gene total count per cell
T = load(['data/Simulated_Branched_Random_Walk_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
N_cell = size(T,2);
N_per_gene = sum(T{:,:},2);
idx_g = find(N_per_gene>N_cell);

% Comupte nearest neighbors for genes with more than 1 UMI per cell on average
for my_norm = My_norm
    if strcmp(My_norm{:},'SanityErrorbar')
% For Sanity distance with errorbar, run:
% Sanity_distance -f data/Simulated_Branched_Random_Walk_UMI_counts.txt -s2n 0 -err 1
% with cut off for genes with more that 1 count per cell on average
        my_D = load('data/cell_cell_distance_with_errorbar.txt');
    else
        load(['data/Simulated_Branched_Random_Walk_' my_norm{:} '_normalization.mat']);
        my_D = pdist(M(idx_g,:)');
    end
    [~,my_knn] = sort(squareform(my_D));
    knn.(my_norm{:}) = my_knn(2:end,:);
end

% Compute fraction of nearest neighbors from distance computed with genes with more than 1 UMI per cell on average
f_nn_mean = zeros(length(K),length(My_norm));
for jj = 1:(length(My_norm));
	for ii = 1:length(K)
		f_nn_mean(ii,jj) = 0;
		for c = 1:N_cell
			N_nn = length(intersect(knn.True(1:K(ii),c),knn.(My_norm{jj})(1:K(ii),c)));
			f_nn_mean(ii,jj) = f_nn_mean(ii,jj) + N_nn/K(ii);
		end
		f_nn_mean(ii,jj) = f_nn_mean(ii,jj)/N_cell;
	end
end
	
% Plot fraction of correct nearest neighbors
subplot(1,3,2)
plot(K',f_nn_mean,'.-');
grid on;
set(gca,'ColorOrder',my_colors);
axis([0 max(K) 0 1])
xlabel('Number of nearest neighbors')
set(gca,'yticklabel',[])
title('Genes with at least 1 UMI per cell','FontWeight','normal')


% Plot Legend :
leg = My_norm;
leg{6} = 'Sanity (without error bars)';
leg{11} = 'Sanity (with error bars)';

y = [6:-1:1 6:-1:2];
x = [1*ones(1,6) 4*ones(1,5)];
my_colors = my_colors(1:length(My_norm),:);

subplot(1,3,3)
scatter(x,y,100,my_colors,'s','filled')
xlim([1 6])
ylim([0.5 6.5])
hold on
for ii = 1:length(y)
	text(x(ii)+.3,y(ii),leg{ii},'FontSize',10)
end
axis off

dim = [24 6];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'figure_7','-dpdf');
