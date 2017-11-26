function plot_ks(spikes_binned, lambdaEsts)

figure(); clf; hold on;
h = [];

for i=1:numel(lambdaEsts)
    
    lambdaEst = lambdaEsts{i};
    timestep = 1;
    lambdaInt = 0;
    j = 0;
    % KS = [];

    for t=1:length(spikes_binned)
        lambdaInt = lambdaInt + lambdaEst(t)*timestep;
        if (spikes_binned(t))
            j = j + 1;
            KS(j) = 1-exp(-lambdaInt);
            lambdaInt = 0;
        end
    end

    KSSorted = sort(KS);
    N = length(KSSorted);

    % plot KS plots
    h(i) = plot(([1:N]-.5)/N, KSSorted,...
            'DisplayName', ['Model' num2str(i)]);

    ks_stat = max(abs(KSSorted - ([1:N]-.5)/N))
    
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