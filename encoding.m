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

