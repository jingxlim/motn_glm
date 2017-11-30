%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutimodal.m
% 
% This function returns covariates for a multimodal place cell using a
% quadratic GLM (eq 2.4 from Agarwal et al.). History dependence is not
% considered in this function.
% 
% Inputs:
%   hist - how far back in time to include as covariates
%   spikes - vector of spikes for this neuron
%   xN - array of normalized x positions
%   yN - array of normalized y positions
% 
% Outputs:
%   cov - array of covariates
%   aic - Akaike Information Criterion
%   ks - KS value
% 
% Function by: Seony ^.^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cov,aic,ks] = multimodal(hist,spikes,xN,yN)
%% covariates
cov = [xN yN xN.^2 yN.^2 xN.*yN];
cov_matrix = hist_dependence(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN);

%% glm
spikes_hist = spikes(hist+1:end);
[b,dev,stats] = glmfit(cov_matrix,spikes_hist,'poisson');

%% analysis
% AIC: Akaike Information Criterion
aic = dev+2*length(b);

% K-S
cov_matrix = [ones(size(cov_matrix,1),1) cov_matrix];
% lambdaEst = exp(sum(b.*cov_matrix));
lambdaEst = gen_lambda(b,cov_matrix);

plot_ks(spikes,lambdaEst);

% ^ I think this needs to be in a loop

% timestep = 1;
% lambdaInt = 0;
% j=0;
% for t=1:length(spikes_binned)
%     lambdaInt = lambdaInt + lambdaEst(t)*timestep;
%     if (spikes_binned(t))
%         j = j + 1;
%         ks(j) = 1-exp(-lambdaInt);
%         lambdaInt = 0;
%     end;
% end;
% KSSorted = sort(ks);
% N = length(KSSorted);
% 
% figure('Name','Models of the Neuron: HW6Q4');
% hold on
% plot( ([1:N]-.5)/N, KSSorted, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
% axis( [0 1 0 1] );
% xlabel('Uniform CDF');
% ylabel('Empirical CDF of Rescaled ISIs');
% title('KS Plot with 95% Confidence Intervals');
% set(gca,'fontsize',14)
% 
% ks_stat = max(abs(KSSorted - ([1:N]-.5)/N));

end