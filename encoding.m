%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 2
% Group: Jingxuan Lim, Simon Orozco, Seony Han
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data & generate variables
load('train.mat')
[ISI,XLocAtSpikes,YLocAtSpikes] = generate_new_variables(x,y,spikes);
% generate velocities?
% generate directions?

%%%% notes %%%%
% more variables could be useful but also means more to present

%% plotting raw data
% plot_spiking_times
plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes);
% plot_spiking_velocities maybe?
% plot_spiking_directions maybe?

%% classifying cells
% place
% grid
% multimodal (?)

%% encoding data
% encode_place
% encode_grid
% encode_mm
% ^ idk if those are functions. read below

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
