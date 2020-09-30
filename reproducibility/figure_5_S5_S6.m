clear all; close all;

addpath('scripts')

% Compare the gene pairwise correlation of two methods (figure 5 Panels A,B, and figure S5)

% Get gene pairwise correlation
load(['data/Baron_Sanity_normalization.mat'])
rho{1} = corr(M');
c{1} = rho{1} - 3*tril(ones(size(rho{1})));
c{1} = c{1}(c{1}>=-1);

for norm_2 = {'TPM' 'scImpute' 'RawCounts' 'Deconvolution' 'sctransform'}
	% Load gene expression and compute gene pairwise correlation
	load(['data/Baron_' norm_2{:} '_normalization.mat'])
	rho{2} = corr(M');
	c{2} = rho{2} - 3*tril(ones(size(rho{1})));
	c{2} = c{2}(c{2}>=-1);

	% Get 2d histogram
	Edges = [-1.005:0.01:1.005]';
	[H X] = hist3([c{2} c{1}],'Edges',{Edges Edges});

	% plot
	figure('visible','off')
	imagesc(X{2},X{1},log10(H),'AlphaData',~(H==0))
	my_map = cbrewer('seq','Blues',256);
	colormap(my_map(64:end,:))
	axis xy;
	hcb = colorbar;
	my_ticks = 0:floor( max(max(log10(H))) );
	for i = 1:length(my_ticks)
		my_ticklabels{i} = ['10^' num2str(my_ticks(i))];
	end
	set(hcb,'Ticks',my_ticks,'TickLabels',my_ticklabels);
	set(get(hcb,'Title') ,'String','nr. of gene pairs');
	hold on;
	plot([-1 0; 1 0],[0 -1; 0 1],'r:')
	axis([-1 1 -1 1])
	xlabel(['Sanity gene pairwise correlation'])
	ylabel([norm_2{:} ' gene pairwise correlation'])

	% print
	dim = [12 9];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
	if any(strcmp(norm_2{:},{'TPM' 'scImpute'}))
		print(gcf,['Fig/figure_5_Sanity_' norm_2{:}],'-dpdf');
	else
		print(gcf,['Fig/figure_S5_Sanity_' norm_2{:}],'-dpdf');
	end
end

% Compare the gene pairwise correlation for high disagreements (figure 5 Panels C-E and figure S6)

% Get raw UMI Counts
T = readtable(['data/Baron_UMI_counts.txt'],'ReadRowNames',1,'delimiter', '\t');
UMI_counts = T{:,:};
[N_gene,N_cell] = size(UMI_counts);

for norm_2 = {'MAGIC' 'DCA' 'SAVER' 'scVI'}
	switch norm_2{:}
	case 'MAGIC'
		th_norm_2{1} = [.975 1.0];
		th_norm_2{2} = [-.3 .3];
		th_norm_1{1} = [-.03 .005];
		th_norm_1{2} = [.6 .93];
	case 'DCA'
		th_conc{1} = [.97 1.0];
		th_conc{2} = [0 .3];
		th_sanity{1} = [-0.05 .005];
		th_sanity{2} = [.55 .75];
	case 'SAVER'
		th_conc{1} = [.88 .97];
		th_conc{2} = [0 .3];
		th_sanity{1} = [-.05 .005];
		th_sanity{2} = [.7 .89];
	case 'scVI'
		th_conc{1} = [.935 1.0];
		th_conc{2} = [-.4 .3];
		th_sanity{1} = [-.03 .005];
		th_sanity{2} = [.6 .93];
	end

	% Load gene expression and compute gene pairwise correlation
	load(['data/Baron_' norm_2{:} '_normalization.mat'])
	rho{2} = corr(M');
	c{2} = rho{2} - 3*tril(ones(size(rho{1})));
	c{2} = c{2}(c{2}>=-1);

	% Get 2d histogram
	Edges = [-1.005:0.01:1.005]';
	[H X] = hist3([c{2} c{1}],'Edges',{Edges Edges});

	% Plot
	figure('visible','off')
	imagesc(X{2},X{1},log10(H),'AlphaData',~(H==0))
	my_map = cbrewer('seq','Blues',256); 
	colormap(my_map(64:end,:))
	axis xy;
	hcb = colorbar;
	my_ticks = 0:floor( max(max(log10(H))) );
	for i = 1:length(my_ticks)
		my_ticklabels{i} = ['10^' num2str(my_ticks(i))];
	end
	set(hcb,'Ticks',my_ticks,'TickLabels',my_ticklabels);
	set(get(hcb,'Title') ,'String','nb. of gene pairs');
	axis([-1 1 -1 1])
	hold on;
	plot(th_norm_1{1}([1 2 2 1 1]),th_norm_2{1}([1 1 2 2 1]),'r-');
	plot(th_norm_1{2}([1 2 2 1 1]),th_norm_2{2}([1 1 2 2 1]),'m-');
	xlabel(['Sanity gene pairwise correlation'])
	ylabel([norm_2{:} ' gene pairwise correlation'])

	% Print
	dim = [12 9];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
	if strcmp(norm_2{:},'MAGIC')
		print(gcf,['Fig/figure_5_Sanity_' norm_2{:}],'-dpdf');
	else
		print(gcf,['Fig/figure_S6_Sanity_' norm_2{:}],'-dpdf');
	end

	% Look at gene pair with high disagrament between the 2 normlisation methods
	for k=1:2
		close all;
		my_map = [];

		% Get index satisfying thresholds
		[i,j] = find( rho{2}>th_norm_2{k}(1) & rho{2}<th_norm_2{k}(2) & rho{1}>th_norm_1{k}(1) & rho{1}<th_norm_1{k}(2) );

		% Remove duplicates
		tmp = sort([i j]')';
		[i,idx]  = sort(tmp(:,1));
		j = tmp(idx,2);
		tmp = unique([i j],'rows');
		i = tmp(:,1);
		j = tmp(:,2);

		% sort pair of genes s.t. count i < count j
		for g = 1:length(i)
			if max(UMI_counts(i(g),:)) > max(UMI_counts(j(g),:))
				tmp = i(g);
				i(g) = j(g);
				j(g) = tmp;
			end
		end

		% Fill 2d-hist;
		H_count = [];
		for g = 1:length(i)
			for l = 1:N_cell
				c_i = UMI_counts(i(g),l)+1;
				c_j = UMI_counts(j(g),l)+1;
				if c_i > size(H_count,1) || c_j > size(H_count,2)
					H_count(c_i,c_j) = 1;
				else
					H_count(c_i,c_j) =  H_count(c_i,c_j) + 1;
				end
			end
		end

		% Sum all count in bin >= 10 in bin 10
		H_count(11,:) = sum(H_count(11:end,:),1);
		H_count(12:end,:) = [];

		% Sum all count in bin >= 30 in bin 30
		H_count(:,31) = sum(H_count(:,31:end),2);
		H_count(:,32:end) = [];

		% Plot
		figure('visible','off')
		imagesc(0:(size(H_count,2)-1),0:(size(H_count,1)-1),log10(H_count),'AlphaData',~(H_count==0));
		axis xy;
		if k==2
			my_map = cbrewer('seq','Purples',256);
			my_map(:,1) = my_map(:,3);

			tmp = cbrewer('div','PiYG',512);
			my_map = flipud( tmp(1:256,:) );
		else
			my_map = cbrewer('seq','Reds',256);
		end
		colormap(my_map(64:end,:))
		hcb = colorbar;
		my_ticks = 0:2:floor( max(max(log10(H_count))) );
		for t = 1:length(my_ticks) 
			my_ticklabels{t} = ['10^' num2str(my_ticks(t))];
		end
		set(hcb,'Ticks',my_ticks,'TickLabels',my_ticklabels);
		set(get(hcb,'Title') ,'String','nb. of cells');
		xtick = 0:5:30;
		xticklabels = {'0','5','10','15' '20','25','30'};
		if max(max(UMI_counts(unique(j),:))) > 30
			xticklabels{end} = '30+';
		end
		ytick = 0:5:10;
		yticklabels = {'0','5','10'};	
		if max(max(UMI_counts(unique(i),:))) > 10
			yticklabels{end} = '10+';
		end
		set(gca,'XTick',xtick,'YTick',ytick,'XTickLabel',xticklabels,'YTickLabel',yticklabels);
		xlim([-.5 30.5])
		ylim([-.5 10.5]);
		xlabel('Count gene i')
		ylabel('Count gene j')
		if k==1
			title(['\color{red}' num2str(length(i)) ' pair of genes'])
		else
			title(['\color{magenta}' num2str(length(i)) ' pair of genes'])
		end

		% Print
		dim = [12 4.5];
		set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
		if strcmp(norm_2{:},'MAGIC')
			print(gcf,['Fig/figure_5_count_hist_Sanity_' norm_2{:} '_' num2str(k)],'-dpdf');
		else
			print(gcf,['Fig/figure_S6_count_hist_Sanity_' norm_2{:} '_' num2str(k)],'-dpdf');
		end
	end
end

