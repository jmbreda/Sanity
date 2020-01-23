clear all; close all;

addpath('scripts')

% define colors
my_colors = flipud(cbrewer('div','RdBu',513));
my_grey = .2*ones(1,3);


% Compare the gene pairwise correlation of two methods (Panels A,B)	
norm_1 = 'Sanity';

% Get gene pairwise correlation
load(['data/Baron_' norm_1 '_normalization.mat'])
rho{1} = corr(M');
c{1} = rho{1} - 3*tril(ones(size(rho{1})));
c{1} = c{1}(c{1}>=-1);

for norm_2 = {'TPM' 'scVI'}
	% Get gene pairwise correlation
	load(['data/Baron_' norm_2{:} '_normalization.mat'])
	rho{2} = corr(M');
	c{2} = rho{2} - 3*tril(ones(size(rho{1})));
	c{2} = c{2}(c{2}>=-1);

	% Get 3d histogram
	Edges = [-1.005:0.01:1.005]';
	[H X] = hist3([c{2} c{1}],'Edges',{Edges Edges});

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

	hold on;
	plot([-1 0; 1 0],[0 -1; 0 1],'r:')

	axis([-1 1 -1 1])
	xlabel([norm_1 ' gene pairwise correlation'])
	ylabel([norm_2{:} ' gene pairwise correlation'])

	dim = [12 9];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
	print(gcf,['Fig/figure_4_' norm_1 '_' norm_2{:}],'-dpdf');
end




% Compare the gene pairwise correlation for high disagreements (Panels C-E)

% Get UMI Counts
T = readtable(['data/Baron_UMI_counts.txt'],'ReadRowNames',1,'delimiter', '\t');
UMI_counts = T{:,:};
[N_gene,N_cell] = size(UMI_counts);

norm_1 = 'Sanity';

% Get gene pairwise correlation
load(['data/Baron_' norm_1 '_normalization.mat'])
rho{1} = corr(M');
clear c;
c{1} = rho{1} - 3*tril(ones(size(rho{1})));
c{1} = c{1}(c{1}>=-1);

for norm_2 = {'MAGIC'}
	switch norm_2{:}
	case 'MAGIC'
		th_norm_2{1} = [.975 1.0];
		th_norm_2{2} = [-.3 .3];
		th_norm_1{1} = [-.03 .005];
		th_norm_1{2} = [.6 .93];
	end

	% Get gene pairwise correlation
	load(['data/Baron_' norm_2{:} '_normalization.mat'])
	rho{2} = corr(M');
	c{2} = rho{2} - 3*tril(ones(size(rho{1})));
	c{2} = c{2}(c{2}>=-1);

	% Get 3d histogram
	Edges = [-1.005:0.01:1.005]';
	[H X] = hist3([c{2} c{1}],'Edges',{Edges Edges});

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
	xlabel([norm_1 ' gene pairwise correlation'])
	ylabel([norm_2{:} ' gene pairwise correlation'])
	dim = [12 9];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
	print(gcf,['Fig/figure_4_' norm_1 '_' norm_2{:}],'-dpdf');

	for k=1:2
		close all;
		my_map = [];

		% Get index satisfying 
		[i,j] = find( rho{2}>th_norm_2{k}(1) & rho{2}<th_norm_2{k}(2) & rho{1}>th_norm_1{k}(1) & rho{1}<th_norm_1{k}(2) );

		% Remove duplicates
		tmp = sort([i j]')';
		[i,idx]  = sort(tmp(:,1));
		j = tmp(idx,2);
		tmp = unique([i j],'rows');
		i = tmp(:,1);
		j = tmp(:,2);

		% sort pair of genes : first lower max count (line)
		for g = 1:length(i)
			if max(UMI_counts(i(g),:)) > max(UMI_counts(j(g),:))
				tmp = i(g);
				i(g) = j(g);
				j(g) = tmp;
			end
		end

		% Fill 3d-hist;
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
		
		H_count(11,:) = sum(H_count(11:end,:),1);
		H_count(12:end,:) = [];

		H_count(:,31) = sum(H_count(:,31:end),2);
		H_count(:,32:end) = [];

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

		dim = [12 4.5];
		set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
		print(gcf,['Fig/figure_4_count_hist_' norm_1 '_' norm_2{:} '_' num2str(k)],'-dpdf');
	end
end

