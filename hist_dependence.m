%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hist_depenence.m
%
% This function creates the matrix of covariates to be used in glmfit when
% spike history is included in the model.
%
% Inputs: num_backsteps - how far back in time to include as covariates
%         spikes - vector of spikes for this neuron
%         varargin - vectors of the non-spike history covariates to be 
%         included
%
% Output: covariate_matrix - the matrix of covariates. Use this with glmfit
%
% Function by: Simon Orozco
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function covariate_matrix = hist_dependence(num_backsteps,spikes,varargin)

n_cov = length(varargin); % Number of non-spike covariates
tot_cov = n_cov + num_backsteps; % Total number of covariates

% Preallocate for speed
covariate_matrix = NaN(length(varargin{1})-num_backsteps,tot_cov);

% Add non-spike covariates to matrix
for i = 1:n_cov
    covariate_matrix(:,i) = varargin{i}(num_backsteps+1:end);
end

% Add spike history dependent covariates to matrix
for n = 1:num_backsteps
    covariate_matrix(:,i+n) = spikes((num_backsteps+1-(n-1)):end-(n-1));
end

end

