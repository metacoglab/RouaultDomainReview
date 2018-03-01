% Demonstration of hierarchical model fits
%
% SF 2014

clear all
HMMpath = '~/Documents/HMM';
addpath(HMMpath);
addpath('~/Dropbox/Utils/meta_d/')
addpath('~/Dropbox/Utils/graphics/export_fig/')
savePlots = 0;
simulate = 0;

Ntrials = [50 400];
Nsub = 100;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;
group_mratio = 0.8;
type1_sigma = 0.2;
rho = 0.5;
type2_sigma = 0.5;

h1 = figure;
set(gcf, 'Position', [200 200 800 600])

for nsim = 1:2
    
    if simulate
        
        N = Ntrials(nsim);
        for i = 1:Nsub
            
            % Generate Mratios for this subject
            bigSigma = [type2_sigma^2 rho.*type2_sigma^2; rho.*type2_sigma^2 type2_sigma^2];
            mratios(i,:) = mvnrnd([group_mratio group_mratio], bigSigma);
            
            %% Task 1
            % Generate dprime
            d = normrnd(group_d, type1_sigma);
            metad = mratios(i,1).*d;
            
            % Generate data
            sim = metad_sim(d, metad, c, c1, c2, N);
            
            nR_S1(1).counts{i} = sim.nR_S1;
            nR_S2(1).counts{i} = sim.nR_S2;
            
            % MLE point estimate fit with padding
            mlefit = fit_meta_d_MLE(sim.nR_S1 + 1/8, sim.nR_S2 + 1/8);
            MLE(1).mratio(i) = mlefit.M_ratio;
            
            %% Task 2
            d = normrnd(group_d, type1_sigma);
            metad = mratios(i,2).*d;
            
            % Generate data
            sim = metad_sim(d, metad, c, c1, c2, N);
            
            nR_S1(2).counts{i} = sim.nR_S1;
            nR_S2(2).counts{i} = sim.nR_S2;
            
            % MLE point estimate fit with padding
            mlefit = fit_meta_d_MLE(sim.nR_S1 + 1/8, sim.nR_S2 + 1/8);
            MLE(2).mratio(i) = mlefit.M_ratio;
        end
        
        % Fit group data all at once
        mcmc_params = fit_meta_d_params;
        fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
        plotSamples(fit.mcmc.samples.mu_logMratio(:,:,1))
        plotSamples(fit.mcmc.samples.mu_logMratio(:,:,2))
        
    else
        
        load(['corrsim' num2str(nsim) '.mat']);
        
    end
    
    figure(h1)
    subplot(2,2,1+((nsim-1)*2))
    plot(MLE(1).mratio, MLE(2).mratio, 'o ', 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 2)
    [rval pval] = corrcoef(MLE(1).mratio, MLE(2).mratio);
    classical_rval(nsim) = rval(2,1);
    classical_pval(nsim) = pval(2,1);
    lsline
    set(gca, 'FontSize', 14, 'XLim', [-0.2 2.5], 'YLim', [-0.2 2.5])
    xlabel('meta-d''/d'' task A')
    ylabel('meta-d''/d'' task B')
    text(1.5,0.2, sprintf('r = %0.2g', classical_rval(nsim)), 'FontSize', 14)
    text(1.5,0,sprintf('p = %0.3g', classical_pval(nsim)), 'FontSize', 14)
    box off
    
    subplot(2,2,2+((nsim-1)*2))
    h= histogram(fit.mcmc.samples.rho(:), 'Normalization', 'probability');
    xlabel('\rho');
    ylabel('Posterior density');
    line([rho rho],[0 max(h.Values)+0.015], 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
    ci = calc_CI(fit.mcmc.samples.rho(:));
    line([ci(1) ci(2)],[0.002 0.002], 'LineWidth', 3, 'Color', [1 1 1])
    box off
    set(gca, 'FontSize', 14, 'XLim', [-1 1])
    
    save(['corrsim' num2str(nsim) '.mat'], 'fit', 'MLE');
    
end

if savePlots
    export_fig('simulatedMetaCorrelation.pdf', '-pdf', '-transparent', '-painters', '-nocrop', h1)
end
