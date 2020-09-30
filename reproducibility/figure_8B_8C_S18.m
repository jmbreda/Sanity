clear all; close all; clc;

addpath('scripts')
addpath('scripts/nmi')

My_clustering = {'ward' 'kmeans' 'louvain'};
Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno_Embryo' 'LaManno_ES' 'LaManno_MouseEmbryo'};
My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};
Measure = {'MI','ARI','FM','Jaccard'};

% define colors
my_colors = load('data/my_colors.txt');

for a = 1:length(My_clustering)
	for d = 1:length(Datasets)
		% Get cell type annotation
		Celltype = textread(['data/' Datasets{d} '_Celltype.txt'],'%s\n');
		% Remove NA: get Clusterings on clustered cells only (i.e. exclude cells annotated NA)
		idx_na = strcmp(Celltype,'NA');
		Celltype = Celltype(~idx_na);
		[~,~,G] = unique(Celltype);

		for n = 1:length(My_norm);

			% load normalised gene expression matrix
			load(['data/' Datasets{d} '_' My_norm{n} '_normalization.mat']);

			switch My_clustering{a}
			case 'ward'
				% Get Euclidean distance
				D = pdist(M(:,~idx_na)');

				% Get tree
				tree = linkage(D,'ward');

				% Get clusters
				K = cluster(tree,'maxclust',max(G));
			case 'kmeans'
				% get kmeans clusters
				K = kmeans(M(:,~idx_na)',max(G),'Replicates',100);
			case 'louvain'
				% Get Euclidean distance
				D = pdist(M(:,~idx_na)');

				% Compute nearest neighbors
				[~,knn] = sort(squareform(D));
				knn = knn(2:end,:);
				N_cell = size(knn,2);

				% print 30 nearest neighbors of each cell for Louvain
				fid = fopen(['data/' Datasets{d} '_' my_norm{:} '_knn30.txt'],'w');
				for c = 1:N_cell
					for k=1:K
						fprintf(fid,'%u %u\n',c-1,knn(k,c)-1);
					end
				end
				fclose(fid);
				% run Louvain clustering algorithm of that file
				% Version 0.2 Based on the article "Fast unfolding of community hierarchies in large networks" Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
				% https://sites.google.com/site/findcommunities/
				% From Louvain_community make data/[Datasets{d}]_[my_norm]_Louvain_tree_knn30.txt
				% From Louvain_hierarchy make data/[Datasets{d}]_[my_norm]_Louvain_hierarchy_knn30.txt
				% Read clusters from tree
				tree_file = ['data/' Datasets{d} '_' my_norm{:} '_Louvain_tree_knn30.txt'];
				hierarchy_file = ['data/' Datasets{d} '_' my_norm{:} '_Louvain_hierarchy_knn30.txt'];
                my_tree = load(tree_file);
                [stat,out] = system(['cut -f3 -d" " ' hierarchy_file '|tail -n+2']);
                N = str2num(out);
                K = [];
                K(:,1) = my_tree(1:N(1),2);
                for ii = 2:length(N-1)
                    idx_tree = ( sum(N(1:ii-1)) + 1 ):( sum(N(1:ii-1)) + N(ii) );
                    my_keys = my_tree(idx_tree,1);
                    my_vals = my_tree(idx_tree,2);
                    if ~all(my_keys==my_vals)
                        K(:,ii) = -1*ones(N(1),1);
                        for jj = 1:length(my_keys)
                            K(K(:,ii-1)==my_keys(jj),ii) = my_vals(jj);
                        end
                    end
                end
                K = K(:,end)+1;
			end
			
			% Get normalized mutual information, ARI, FM, Jaccard
			[mi,ari,fm,jaccard] = get_mi_ari_fm_jaccard(G,K);

			% save in score 
			score(a,n,d,1) = mi;
			score(a,n,d,2) = ari;
			score(a,n,d,3) = fm;
			score(a,n,d,4) = jaccard;
		end
	end
end

% Get number of time each normalization method obtains the best score
n_best = zeros(1,length(My_norm));
for m = 1:length(Measure)
	for d = 1:length(Datasets)
		for a = 1:length(My_clustering)
			n_best(1,:) = n_best(1,:) + double(score(a,:,d,m)==max(score(a,:,d,m)));
		end
	end
end

% plot n_best as bar
figure('visible','off');
hold on;
for n = 1:length(My_norm)
	bar(n,n_best(n),'FaceColor',my_colors(n,:));
end
box on
axis([0.5 length(My_norm)+.5 0 45])
set(gca,'XTick',1:length(My_norm),'XTickLabel',My_norm,'XTickLabelRotation',60)
ylabel('N')
title('Number of times each method performs best','FontWeight','normal')

dim = [10 8];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'figure_8B','-dpdf');

% Compute ratio of each score to the best score in each combination
ratio_to_best = -1*ones(length(Measure)*length(Datasets)*length(My_clustering),length(Measure));
k = 0;
for m = 1:length(Measure)
	for d = 1:length(Datasets)
		for a = 1:length(My_clustering)
			k=k+1;
			for n = 1:length(My_norm)
				ratio_to_best(k,n) = score(a,n,d,m)/max(score(a,:,d,m));
			end
		end
	end
end

% get quantiles on ratio_to_best
q = [.05 .25 .5 .75 .95];
for i = 1:5
	Quant(i,:) = quantile(ratio_to_best,q(i));
end

% plot ratio_to_best as box plot
figure('visible','off');
hold on;
for n = 1:length(My_norm)
	% 5% - 95%
	plot([n;n],Quant([1 5],n),'linewidth',1,'color',my_colors(n,:))
	% 25% - 75%
	plot([n;n],Quant([2 4],n),'linewidth',8,'color',my_colors(n,:))
	% 50%
	plot(n,Quant(3,n),'o','color',my_colors(n,:),'MarkerFaceColor',[1 1 1])
	plot(n,Quant(3,n),'.','color',my_colors(n,:))
end
axis([0.5 length(My_norm)+.5 0 1])
set(gca,'XTick',1:length(My_norm),'XTickLabel',My_norm,'XTickLabelRotation',60)
ylabel('score/score best method')
title('Dist. performance relative to best performance','FontWeight','normal');
box on

dim = [10 8];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'figure_8C','-dpdf');

% get sorting index based on mean score
[~,idx_n] = sort( mean(mean(mean( score(:,:,:,:),4),3),1), 'descend');

% plot all scores as bar
figure('visible','off'); f=1;
for d = 1:length(Datasets)
	for a = 1:length(My_clustering)
		my_score = reshape(score(a,idx_n,d,:),length(My_norm),length(Measure));

		subplot(length(Datasets)+1,length(My_clustering),f);f=f+1;
		bar(1:length(Measure),my_score,1,'EdgeColor','none')
		set(gca,'ColorOrder',my_colors(idx_n,:));
		box off
		axis([.5 length(Measure)+.5 0 1])
		set(gca,'xtick',[])
		if d==1
			title(Algo_label{a},'FontWeight','normal')
		end
		if a==1
			ylabel(Datasets_label{d})
		end
		if d==length(Datasets)
			set(gca,'xtick',1:length(Measure),'xticklabel',Measures,'xticklabelrotation',0)
		end
	end
end

% legend
subplot(length(Datasets)+1,1,length(Datasets)+1)
y = [2*ones(length(My_norm)/2,1);ones(length(My_norm)/2,1)];
x = repmat([1:length(My_norm)/2]',2,1);
scatter(x,y,100,my_colors(idx_n,:),'o','filled')
hold on
for ii = 1:length(y)
	text(x(ii)+.1,y(ii),My_norm{idx_n(ii)},'FontSize',10)
end
axis([1 length(My_norm)/2+.5 1 2.5])
axis off

dim = [6*length(My_clustering) 3*length(Datasets)];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'figure_S18','-dpdf');
