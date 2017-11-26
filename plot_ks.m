% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot_ks.m
% -------------------------------------------------------------------------
% plot_ks takes in the spikes of a single neuron (spikes_binned) and
% several models (lambdaEsts) and plot the Kolmogorov-Smirnov plots and
% computes the KS static for each of them.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function plot_ks(spikess, lambdaEsts)

figure(); clf; hold on;
h = [];

for i=1:numel(lambdaEsts)
    disp(['Model: ' num2str(i)])
    
    lambdaEst = lambdaEsts{i};
    spikes = spikess{i};
    timestep = 1;
    lambdaInt = 0;
    j = 0;
    % KS = [];

    for t=1:length(spikes)
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes(t))
            j = j + 1;
            KS(j) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted = sort(KS);
    N = length(KSSorted)

    % plot KS plots
    ks_stat = max(abs(KSSorted - ([1:N]-.5)/N));
    h(i) = plot(([1:N]-.5)/N, KSSorted,...
            'DisplayName', ['KS = ' num2str(ks_stat)]);
    
end

% plot "perfect" models
plot(0:.01:1,0:.01:1, 'k',...
    0:.01:1, [0:.01:1]+1.36/sqrt(N),'k--',...
    0:.01:1,[0:.01:1]-1.36/sqrt(N), 'k--')


axis( [0 1 0 1] );
xlabel('Uniform CDF');
ylabel('Empirical CDF of Rescaled ISIs');
title('KS Plot with 95% Confidence Intervals');
legend(h, 'Location', 'northwest');