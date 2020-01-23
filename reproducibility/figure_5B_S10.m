clear all; close all; clc;

addpath('scripts')
addpath('scripts/nmi')

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo'};
My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};

% define colors
my_colors = cbrewer('qual','Set1',length(My_norm));
tmp_colors = cbrewer('qual','Dark2',8);
my_colors(6,:) = mean([my_colors(6,:); tmp_colors(6,:)]);

my_label = [Datasets 'Average'];

fig_name.ward.NMI =   '5B';
fig_name.ward.ARI =   'S10A';
fig_name.kmeans.NMI = 'S10B';
fig_name.kmeans.ARI = 'S10C';

for my_clustering = {'ward' 'kmeans'}

	Clust.NMI = nan(length(Datasets),length(My_norm));
	Clust.ARI = nan(length(Datasets),length(My_norm));

	for k = 1:length(Datasets)
		dataset = Datasets{k};

		% Get cell type annotation
		Celltype = textread(['data/' dataset '_Celltype.txt'],'%s\n');
		% remove NA
		idx_na = strcmp(Celltype,'NA');
		Celltype = Celltype(~idx_na);
		[~,~,G] = unique(Celltype);

		% Get Clusterings on clustered cells only (i.e. exclude not NA)
		for n=1:length(My_norm);
			load(['data/' dataset '_' My_norm{n} '_normalization.mat'])

			if strcmp(my_clustering{:},'ward')
				% Get Euclidean distance
				D = pdist(M(:,~idx_na)');
				% Get tree
				tree = linkage(D,'ward');
				% Get clusters
				K = cluster(tree,'maxclust',max(G));
			else
				K = kmeans(M(:,~idx_na)',max(G),'Replicates',100);
			end
			
			% Get normalized mutual information
			Clust.NMI(k,n) = nmi(G,K);

			% Get Confision Matrix " cm(1,1) = TP, cm(1,2) = FN, cm(2,1) = FP, cm(2,2) = TN
			cm = zeros(2,2);
			N_cell = length(K);
			for i = 1:(N_cell-1)
				for j = (i+1):N_cell
					Truth = G(i)==G(j);
					Seen  = K(i)==K(j);
					cm(-Truth+2,-Seen+2) = cm(-Truth+2,-Seen+2)+1;
				end
			end

			% Get adjusted Rand index
			Clust.ARI(k,n) = 2*(cm(1,1)*cm(2,2) - cm(1,2)*cm(2,1))/( (cm(1,1)+cm(1,2))*(cm(1,2)+cm(2,2)) + (cm(1,1)+cm(2,1))*(cm(2,1)+cm(2,2)) );

		end
	end


	% Bar Plots
	my_grey = [.2 .2 .2];
	for c = fieldnames(Clust)'

		% Get averaged score over datasets for ranking
		Clust.(c{:})(length(Datasets)+1,:) = mean(Clust.(c{:}),1,'omitnan');
		[~,idx_norm_sorting] = sort(Clust.(c{:})(end,:),'descend');

		figure('visible','off')
		bar(abs(Clust.(c{:})(:,idx_norm_sorting)),1);
		colormap(my_colors(idx_norm_sorting,:));
		legend(My_norm(idx_norm_sorting),'location','EastOutside');	
		axis([0.4 7.6 0 1])
		set(gca,'xticklabel',my_label,'xticklabelrotation',20)
		xlabel('Datasets')
		ylabel('Similarity score')

		dim = [24 10];
		set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
		'PaperPosition',[0 0 dim],'PaperSize',[dim])
		print(gcf,['Fig/figure_' fig_name.(my_clustering{:}).(c{:})],'-dpdf');
	end
end
