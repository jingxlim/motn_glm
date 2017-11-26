%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 2
% Group: Jingxuan Lim, Simon Orozco, Seony Han
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start anew
clear all; % clear previous variables
close all; % close previous plots

%% load data & generate variables

load('train.mat')
[ISI,XLocAtSpikes,YLocAtSpikes,Vx,Vy,dir] = generate_new_variables(xN,yN,spikes_binned);
% generate directions?

ISI_threshold = 550;
%%%% notes %%%%
% more variables could be useful but also means more to present

%% plotting raw data
% plot_spiking_times
plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes,'subplot');
ISIs = plot_ISIs(spikes_binned,ISI_threshold,2);

% plot_spiking_velocities maybe?
% plot_spiking_directions maybe?

%% classifying cells
% This needs to be automated, and it turns out that is a hard problem (at
% least I struggled with it -Simon)

% place
% grid
% multimodal (?) (I think this is where the history dependence could prove
% to be a helpful covariate -Simon)

%% encoding data
% encode_place
% encode_grid
% encode_mm
% ^ idk if those are functions. read below

for i=6:8
    disp(['Working on neuron ' num2str(i) ' ...'])
    spikes = spikes_binned(:,i);
    
    % DEFINE MODELS    
    % Model 1
    [b,dev,stats] = glmfit([xN yN xN.^2 yN.^2],spikes,'poisson');
    lambdaEst{1} = exp(b(1)+b(2)*xN+ b(3)*yN+b(4)*xN.^2+ b(5)*yN.^2);
    spikess{1} = spikes;
    % Model 2
    [b,dev,stats] = glmfit([xN yN xN.^2 yN.^2 xN.*yN],spikes,'poisson');
    lambdaEst{2} = exp(b(1)+b(2)*xN+ b(3)*yN+b(4)*xN.^2+ b(5)*yN.^2+b(6)*xN.*yN);
    spikess{2} = spikes;
    % Model 3
    [b,dev,stats] = glmfit([xN yN],spikes,'poisson');
    lambdaEst{3} = exp(b(1)+b(2)*xN+ b(3)*yN);
    spikess{3} = spikes;
    % Model 4: history dependence
    hist = 10;
    spikes_m4 = spikes(hist+1:end);
    covar_m = hist_dependence(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN);
    [b,dev,stats] = glmfit(covar_m,spikes_m4,'poisson');
    lambdaEst{4} = exp(b(1)+b(2)*covar_m(:,1)+ b(3)*covar_m(:,2)+...
                       b(4)*covar_m(:,3)+b(5)*covar_m(:,4)+...
                       b(6)*covar_m(:,5)+b(7)*covar_m(:,6)+...
                       b(8)*covar_m(:,7)+b(9)*covar_m(:,8)+...
                       b(10)*covar_m(:,9)+b(11)*covar_m(:,10)+...
                       b(12)*covar_m(:,11)+b(13)*covar_m(:,12)+...
                       b(14)*covar_m(:,13)+b(15)*covar_m(:,14)+...
                       b(16)*covar_m(:,15));
    spikess{4} = spikes_m4;
    
    plot_ks(spikess,lambdaEst);
end

%%%% notes %%%%
% what needs to happen
%   - glmfit (diff covariates for each cell group)
%   - lambda calc (could do in same fn as w/glmfit) --> calc_glm_lambda
%   - plot lambda

% options for implementation
%   1. call calc_glm_lambda and plot_model for each cell group
%   2. make separate encoding functions for each cell group. Those
%   functions would call calc_glm_lambda and plot_model

% other things
%   - covariates depending on factors besides position would need different
%   plot_model

% These functions could all call the same glm/lambda fn (calc_glm_lambda)
%   - inputs: covariates, spikes_binned
%   - outputs: b, dev, stats --> double check how many of these we acutally
%   need outside of the glm

%% validate covariate choices

