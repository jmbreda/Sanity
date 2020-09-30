clear all; close all; clc;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','Deconvolution','MAGIC','Sanity','SAVER','scImpute','sctransform','scVI'};
Datasets = {'Zeisel' 'LaManno_Embryo' 'LaManno_MouseEmbryo'};

% Compute sensitivity and Positive predictive value
for d = 1:length(Datasets)
	
	% Load refence DE values and gene names
	my_file = ['data/' Datasets{d} '_Reference_differencial_expression.txt'];
	RefMarker = readtable(my_file,'ReadRowNames',1,'delimiter','\t');

    % Get cell type annotation
    Celltype = textread(['data/' Datasets{d} '_Celltype.txt'],'%s\n');
    [u_Celltype,~,G] = unique(Celltype);

    % put NA cells (unclustered) as cluster 0
    idx_na = strcmp(Celltype,'NA');
    if ~isempty(idx_na)
        k_na = unique(G(idx_na));
        G(idx_na) = 0;
        for k = k_na:(max(G)-1)
            G(G==k+1) = k;
        end
        u_Celltype(k_na) = [];
    end

	for n = 1:length(My_norm)
		load(['data/' d{:} '_' My_norm{n} '_normalization.mat']);

		% Compute zscore per gene per celtype
		Z = zeros(size(M,1),max(G));
		for k = 1:max(G)
			N_k = sum(G==k);
			N_not_k = sum(G~=k);
			Z(:,k) = (mean(M(:,G==k),2) - mean(M(:,G~=k),2))./sqrt(std(M(:,G==k),0,2).^2/N_k +        std(M(:,G~=k),0,2).^2/N_not_k);
		end
	
		% Get gene names
		T = readtable(['data/' Datasets{d} '_UMI_counts.txt'],'ReadRowNames',1,'delimiter','\t');
		GeneID = T.Properties.RowNames;
		clear T;

		% Make Z a table with gene as rows and cell types as colums
		Z = array2table(Z);
		Z.Properties.RowNames = GeneID;
		Z.Properties.VariableNames = strrep(u_Celltype,'-','_');
		
		% Intersect the genes of the reference table with all genes of the Z-score table
		[~,tmp,idx_z] = intersect(RefMarker.Properties.RowNames,Z.Properties.RowNames,'stable');
		if ~all(tmp'==1:size(RefMarker.Properties.RowNames))
			disp('error')
			break;
		end

		% Intersect the celltypes of the reference and the Z-score
		[~,idx_type_ref,idx_type_z] = intersect(RefMarker.Properties.VariableNames,Z.Properties.VariableNames);
		Z = Z{idx_z,idx_type_z};
		Z = Z(:);
		I_Ref = RefMarker{:,idx_type_ref};
		I_Ref = I_Ref(:);
		N = [1:length(I_Ref)]';
		
		% Compute sensitivity and Positive predictive value
		[Z,idx_sort] = sort(Z,'descend');
		sensitivity.(Datasets{d}).(My_norm{n}) = cumsum(I_Ref(idx_sort))./sum(I_Ref);
		PPV.(Datasets{d}).(My_norm{n}) = cumsum(I_Ref(idx_sort))./N;

		% Dot at P_Z = 0.95
		P = .5*(1-erf(-Z/sqrt(2)));
		P = P(:);
		P_th = abs(P-0.95);
		IDX_P_th.(Datasets{d}).(My_norm{n}) = min( find( P_th==min(P_th) ) );
	end
end


% define colors
load('data/my_colors.txt')

figure('visible','off');
for d = 1:length(Datasets)
	
	subplot(2,2,d)
	hold on
	% plot sensitivity against Positive predictive value curve
	for n=1:length(My_norm)
		sens = sensitivity.(Datasets{d}).(My_norm{n});
		ppv = PPV.(Datasets{d}).(My_norm{n});
		tmp = uniquetol([sens ppv],5e-4,'ByRows',1);
		sens = tmp(:,1);
		ppv = tmp(:,2);

		plot(sens,ppv,'linewidth',1,'color',my_colors(n,:));
	end

    % Add dot at P_Z = .95
	for n=1:length(My_norm)
		sens = sensitivity.(Datasets{d}).(My_norm{n});
		ppv = PPV.(Datasets{d}).(My_norm{n});
		idx = IDX_P_th.(Datasets{d}).(My_norm{n});

		plot(sens(idx),ppv(idx),'o','color',my_colors(n,:),'Markerfacecolor',my_colors(n,:),          'Markersize',8);
	end
	axis([0 1 0 1])
	xlabel('Sensitivity')
	ylabel('Pos. predictive value')
	title({'Reference differential expression',strrep(Datasets{d},'_',' ')},'FontWeight','normal')
end

% Add Legend
subplot(2,2,4)
C = my_colors(1:length(My_norm),:);
y = [length(My_norm)/2:-1:1 length(My_norm)/2:-1:1];
x = [zeros(1,length(My_norm)/2) 3*ones(1,length(My_norm)/2)];
scatter(x,y,100,C,'s','filled')
hold on
for ii = 1:length(y)
	text(x(ii)+.2,y(ii),My_norm{ii},'FontSize',10)
end
axis([-.2 5.5 1 6])
axis off

% Print
dim = [20 14];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_S24','-dpdf');
