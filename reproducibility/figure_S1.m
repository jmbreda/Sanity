clear all; close all; clc;

Spieces = {'hg19','mm10'};

figure('visible','off'); f=1;
for s = Spieces
	disp(s{:})

	% load gene expression data, compute TPM. (file on Zenodo)
	T = readtable(['data/fantom5_' s{:} '_expression.txt'],'ReadRowNames',1,'delimiter','\t');
	TPM = T{:,:}./sum(T{:,:})*1e6;

	% Compute mean and variance of TPM
	mu = mean(TPM,2);
	sig2 = var(TPM,0,2);

	% plot variance vs. mean scatter
	subplot(2,2,f);f=f+1;
	scatter(mu,sig2,'.');
	hold on
	x = [1e-2,1e5];
	y = [1e-4,1e10];
	plot(x,y,'linewidth',2)
	set(gca,'xscale','log','yscale','log')
	set(gca,'xtick',10.^[-4:5],'ytick',10.^[-4:2:10])
    ylim(y)
	set(gca,'XGrid','on','XMinorGrid','off')
	set(gca,'YGrid','on','YMinorGrid','off')
	legend('data','$f(x)=x^{2}$','location','SouthEast','interpreter','latex')
	xlabel('mean TPM')
	ylabel('variance TPM')
	title(s{:},'FontWeight','normal')

	% compute mean and variance of log TPM
	mu = mean(log(TPM),2);
	sig2 = var(log(TPM),0,2);

    % plot variance vs. mean scatter
	subplot(2,2,f);f=f+1;
	scatter(mu,sig2,'.');
	set(gca,'xtick',[0:2:10],'ytick',[0:2:10])
	grid on;
	xlabel('mean log TPM')
	ylabel('variance log TPM')
	title(s{:},'FontWeight','normal')
end
dim = [20 15];
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0  dim],'PaperSize',[dim]);
print(gcf,'Fig/figure_S1','-dpdf');

