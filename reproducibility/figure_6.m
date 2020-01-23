clear all; close all; clc;

addpath('scripts')

My_norm = {'RawCounts','TPM','DCA','MAGIC','Sanity','SAVER','scImpute','scVI'};
Datasets = {'Zeisel' 'LaManno_Embryo' 'LaManno_MouseEmbryo'};


% Compute sensitivity and Positive predictive value
for d = 1:length(Datasets)

	my_file = ['data/' Datasets{d} '_Reference_differencial_expression.txt'];
	RefMarker = readtable(my_file,'ReadRowNames',1,'delimiter','\t');

    % Get cell type annotation
    Celltype = textread(['../data/' Datasets{d} '/Celltype.txt'],'%s\n');
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
		load(['data/' d{:} '_' my_norm{:} '_normalization.mat']);

		%Get zscore per gene per celtype
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
		% Make Z a table with gene names and cell types
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
my_colors = cbrewer('qual','Set1',length(My_norm));
tmp_colors = cbrewer('qual','Dark2',8);
my_colors(6,:) = mean([my_colors(6,:); tmp_colors(6,:)]);

% plot sensitivity against Positive predictive value
for d = 1:length(Datasets)
	
	figure('visible','off');
	hold on
	for n=1:length(My_norm)
		sens = sensitivity.(Datasets{d}).(My_norm{n});
		ppv = PPV.(Datasets{d}).(My_norm{n});
		tmp = uniquetol([sens ppv],5e-4,'ByRows',1);
		sens = tmp(:,1);
		ppv = tmp(:,2);

		plot(sens,ppv,'linewidth',1,'color',my_colors(n,:));
	end

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

	dim = [10 7];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
	'PaperPosition',[0 0 dim],'PaperSize',[dim])
	print(gcf,['Fig/PPV_Sensitivity_' Datasets{d} '_Ref'],'-dpdf');
end
