clear all; close all; clc;

% Get UMI count table
tmp = readtable(['data/SimulatedBaron_UMI_counts.txt'],'ReadRowNames',1,'delimiter', '\t');
T.sim = tmp;
clear tmp;
tmp = readtable(['data/Baron_UMI_counts.txt'],'ReadRowNames',1,'delimiter', '\t');
T.baron = tmp;
clear tmp;

% Count UMI per cell and per gene
for t = {'sim','baron'}
	N_c.(t{:}) = sum(T.(t{:}){:,:},1)';
	N_c.(t{:}) = log10(N_c.(t{:}));
	N_g.(t{:}) = sum(T.(t{:}){:,:},2);
	N_g.(t{:}) = log10(N_g.(t{:}));
end

% Define bins
bins_Nc = linspace(min([N_c.sim; N_c.baron]),1.001*max([N_c.sim; N_c.baron]),42);
dx_Nc = bins_Nc(2)-bins_Nc(1);
bins_Ng = linspace(0,1.001*max([N_g.sim; N_g.baron]),42);
dx_Ng = bins_Ng(2)-bins_Ng(1);

% Get distributions of UMI counts
for t = {'sim','baron'}
	p_Nc.(t{:}) = histc(N_c.(t{:}),bins_Nc);
	p_Nc.(t{:})(end) = [];
	p_Nc.(t{:}) = p_Nc.(t{:})/(sum(p_Nc.(t{:})));

	p_Ng.(t{:}) = histc(N_g.(t{:}),bins_Ng);
	p_Ng.(t{:})(end) = [];
	p_Ng.(t{:}) = p_Ng.(t{:})/(sum(p_Ng.(t{:})));
end
bins_Nc = .5*(bins_Nc(1:end-1) + bins_Nc(2:end))';
bins_Ng = .5*(bins_Ng(1:end-1) + bins_Ng(2:end))';

% Compute variance
for t = {'sim','baron'}
	v.(t{:}) = log10( var(T.(t{:}){:,:},0,2) );
end
bins_v = linspace(min([v.sim;v.baron]),1.001*max([v.sim;v.baron]),42);

% Get distribution of variance
for t = {'sim','baron'}
	p_v.(t{:}) = histc(v.(t{:}),bins_v);
	p_v.(t{:})(end) = [];
	p_v.(t{:}) = p_v.(t{:})/(sum(p_v.(t{:})));
end
bins_v = .5*(bins_v(1:end-1) + bins_v(2:end))';

% Plot the distributions
figure('visible','off')
subplot(1,3,1)
plot(bins_Nc,[p_Nc.sim p_Nc.baron],'linewidth',3)
xlim([min(bins_Nc) max(bins_Nc)]);
xlabel('Total mRNA per cell')
ylabel('Frequency')

xtick = 3:4;
xticklabel = {'10^3' '10^4'};
set(gca,'Xtick',xtick,'XtickLabel',xticklabel);

subplot(1,3,2)
plot(bins_Ng,[p_Ng.sim p_Ng.baron],'linewidth',3)
xlim([0 max(bins_Ng)]);
xlabel('Total mRNA per gene')
ylabel('Frequency')

xtick = 0:6;
xticklabel = {'10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6'};
set(gca,'Xtick',xtick,'XtickLabel',xticklabel);

subplot(1,3,3)
plot(bins_v,[p_v.sim p_v.baron],'linewidth',3)
xlim([min(bins_v) max(bins_v)]);
legend('Simulated','Baron')
xlabel('Gene variance')
ylabel('Frequency')

xtick = -2:2:4;
xticklabel = {'10^{-2}' '10^0' '10^2' '10^4'};
set(gca,'Xtick',xtick,'XtickLabel',xticklabel);

dim = [24 6];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
'PaperPosition',[0 0 dim],'PaperSize',[dim])
print(gcf,['Fig/figure_S11'],'-dpdf');
