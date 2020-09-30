clear all; close all;

addpath('scripts')
My_norm = {'True' 'SanityErrorbar' 'scVI' 'DCA' 'SAVER' 'Sanity' 'scImpute' 'RawCounts' 'Deconvolution' 'sctransform' 'MAGIC'}; % Sorted my integral of fiugre 7

for my_norm = My_norm
	% Compute Euclidean distances
	if strcmp(my_norm{:},'True')
		load('data/Simulated_Branched_Random_Walk.mat');
		N_cell = length(neighbors);
		D = squareform(D_true);
	end
		load(['data/Simulated_Branched_Random_Walk_' my_norm{:} '_normalization.mat'])
		D = squareform(pdist(M'));
	end

	% Print distance matrix
	dlmwrite(['data/Simulated_Branched_Random_Walk_' my_norm{:} '_Euclidean_dist.txt'],D,'delimiter','\t');
end

% Run tsne on the distances matrices
[stat,out] = system('python get_tsne.py');

% Define cell colors
my_colors_149 = distinguishable_colors(N_path);
c = 0;
for i = 1:N_cell
	if mod(i-1,13)==0
		c=c+1;
	end
	my_colors(i,:) = my_colors_149(c,:);
end

% Get figure dimention
p0 = 0.05;
dx = 0.05;
dy = 0.04;
width =  (1 -2*p0 -2*dx)/3;
height = (1 -2*p0 -3*dy)/4;
x = repmat([p0 p0+width+dx p0+2*(width+dx)],4,1)';
y = 1 - repmat([p0 p0+height+dy p0+2*(height+dy) p0+3*(height+dy)],3,1)-height;
d = 7;
dim = [d/width d/height];

figure('visible','off');f=1;
for my_norm = My_norm
	% Load tsne visualisation (from get_tsne.py)
	tsne =load(['data/Simulated_Branched_Random_Walk_' my_norm{:} '_tsne.txt']);

	% Plot
	subplot('Position',[x(f) y(f) width height]);f=f+1;
	scatter(tsne(:,1),tsne(:,2),1,my_colors,'.')	
	axis([min(tsne(:,1)) max(tsne(:,1)) min(tsne(:,2)) max(tsne(:,2))])
	axis off
	title(my_norm{:})
end

% Print
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,'Fig/figure_S12','-dpdf');
