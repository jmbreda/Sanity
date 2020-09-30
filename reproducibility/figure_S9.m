clear all; close all; clc;

addpath('~/Matlab_Scripts/cbrewer/')

% Load TPM normalization
load('data/Simulated_Baron_Independent_Genes_TPM_normalization.mat');
[N_gene,N_cell] = size(M);

% Compute svd
M = M-mean(M,2);
[U,S,V] = svd(M','econ');
coeff = V;
score = U*S;

% Number of principal components
N_pc = [N_gene 1000 200 100 50 30 20 10 5];

% Compute TPM correlation matrix
Corr_pc = zeros(N_gene,N_gene,length(N_pc));
Corr_pc(:,:,1) = corr(M');

% Compute correlation matrix in n first principal components
for ii = 2:length(N_pc)
	Corr_pc(:,:,ii) = corr(score(:,1:N_pc(ii))*coeff(:,1:N_pc(ii))');
end

% Compute true correlation matrix
load('data/Simulated_Baron_Independent_Genes.mat');
C = corr(E');
my_idx = find(triu(ones(N_gene),1));
c{1} = C(my_idx);

figure('visible','off'); f=1;
for ii = 1:length(N_pc)

    tmp = Corr_pc(:,:,ii);
    c{2} = tmp(my_idx);

	% Make 2d histogram
    Edges = [-1.005:0.01:1.005]';
    [H X] = hist3([c{2} c{1}],'Edges',{Edges Edges});

	% Plot
    subplot(3,3,f);
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
    title([num2str(N_pc(ii)) ' pcs'],'FontWeight','normal')
    if f==9
        xlabel(['True Pearson correlation'])
    end
    if f==5
        ylabel(['Inferred Pearson correlation'])
    end
end

dim = [20 15];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S9','-dpdf');
