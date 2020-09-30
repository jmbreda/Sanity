clear all; close all; clc;

Datasets = {'Zeisel' 'Baron' 'Chen' 'LaManno/Embryo' 'LaManno/ES' 'LaManno/MouseEmbryo'};

% Load mean log transcription quotient from Sanity
for d = Datasets
	d_name = strrep(d{:},'/','_');
	load(['data/' d{:} '_Sanity_normalization.mat'])
	mu.(d_name) = mean(M,2);
    mu.(d_name) = mu.(d_name) - log(sum(exp(mu.(d_name))));
end

% Compute distribution of mean ltqs
H_mu = [];
x_mu = [];
n_bin = 30;
for d = 1:length(Datasets)
	d_name = strrep(Datasets{d},'/','_');

	% Define bins
    mu_bins = linspace( min(mu.(d_name)), max(mu.(d_name))*1.01,n_bin)';

	% Compute histogram
	H_mu(:,d) = histc(mu.(d_name),mu_bins);
	
	% Compute bin centers (i.e. mean mu in each bin)
	N_gene = length(mu.(d_name));
    for b = 1:n_bin-1
        idx = find(mu.(d_name) >= mu_bins(b) & mu.(d_name) < mu_bins(b+1));
        x_mu(b,d) = mean(mu.(d_name)(idx));
    end

	% Normalize
    dx = mu_bins(2:end)-mu_bins(1:end-1);
    dx = 1;
	H_mu(1:end-1,d) = H_mu(1:end-1,d)./(dx*sum(H_mu(:,d)));
end
% Delete last bin (empty)
H_mu(end,:) = [];

% Plot mean ltq distribution
figure('visible','off')
sz = 1;
loglog(exp(x_mu),H_mu,'o-','linewidth',sz,'MarkerSize',sz);
xlabel('<Gene transcription quotient>')
ylabel('Frequency')
legend(Datasets,'location','SouthWest')
set(gca,'XTick',10.^[-7:-2])

% Print
dim = [10 7];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S22','-dpdf');


% Plot reverse cumulative of mean tq
figure('visible','off')
for d = Datasets
	d_name = strrep(d{:},'/','_');
    loglog(sort(exp(mu.(d_name))),linspace(1,0,length(mu.(d_name))),'-','linewidth',sz,'MarkerSize',sz);
    hold on
end
set(gca,'XTick',10.^[-7:-2],'YTick',10.^[-5:0])
xlabel('<Gene TQ>')
ylabel('Reverse cdf')
legend(Datasets,'location','EastOutside')

% Print
dim = [15 7];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',         'PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/tq_rev_cdf','-dpdf');

% Pool together all mean ltq vectors
MU = [];
for d = fieldnames(mu)
    MU = [MU; mu.(d_name)];
end
% Total UMI per cell
N = [1e3 3*1e3 1e4 3*1e4 1e5];

% Plot cumulative distribution of gene mean UMI per cell for different given total UMI per cell
figure('visible','off')
leg={};
f_above_1 = [];
y = linspace(1,0,length(MU));
colors = lines(5);
subplot(1,2,1)
for n = 1:length(N)
    leg{end+1} = [num2str(N(n))];
    loglog(sort(exp(MU))*N(n),y,'-','linewidth',sz);
    hold on
    f_above_1(n) = sum( exp(MU)*N(n) > 1 )/length(MU);
end
axis([min(sort(exp(MU))*N(1)) max(sort(exp(MU))*N(end)) min(y) max(y)])
set(gca,'XTick',10.^[-4:4],'YTick',10.^[-5:0])
xlabel('<Gene UMI per cell>')
ylabel('Reverse cdf')
legend(leg,'location','SouthWest')

% Plot fraction of genes with more that 1 UMI per cell on average as a function of the average total UMI per cell.
subplot(1,2,2)
semilogx(N,f_above_1,'k-')
hold on
for n = 1:length(N)
    plot(N(n),f_above_1(n),'o','color',colors(n,:),'MarkerFaceColor',colors(n,:),'MarkerSize',5);
end
set(gca,'XTick',N,'XtickLabel',{'1000','3000','10000','30000','100000'})
xlabel('<UMI per cell>')
ylabel({'fraction of genes with' 'more than 1 UMI per cell'})

% Print
dim = [20 7];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',         'PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S23','-dpdf');
