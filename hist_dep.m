%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hist_dep.m
%
% This function creates the matrix of covariates to be used in glmfit when
% spike history is included in the model.
%
% Inputs: timesteps - an array of time steps counted in reverse used for
%                     history dependence (i.e. 1 is the most recent,
%                     followed by 2...)
%            spikes - vector of spikes for this neuron
%          varargin - vectors of the non-spike history covariates to be 
%                     included
%
% Output: covariate_matrix - the matrix of covariates. Use this with glmfit
%
% Function by: Lim Jing Xuan (adapted from hist_dependence.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spikes_trunc,covar_m] = hist_dep(timesteps,spikes,varargin)

n_cov = length(varargin); % Number of non-spike covariates
n_timesteps = length(timesteps);  % Number of spike convariates
tot_cov = n_cov + n_timesteps; % Total number of covariates

spikes_trunc = spikes(max(timesteps)+1:end);

% Preallocate for speed
covar_m = zeros(length(spikes_trunc),tot_cov);

% Add non-spike covariates to matrix
for i = 1:n_cov
    covar_m(:,i) = varargin{i}(max(timesteps)+1:end);
end

% Add spike history dependent covariates to matrix
for n = timesteps
    covar_m(:,n_cov+n) = spikes(max(timesteps)-n+1:end-n);
end

if length(spikes_trunc) ~= size(covar_m,1);
    warning('Dataset size mis-match')
end

end