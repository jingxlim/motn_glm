% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% gen_lambda.m
% -------------------------------------------------------------------------
%
% This function generates the conditional intensity of a model given the
% b coefficients and the covariates with which the coefficients are
% derived from using glmfit.
%
% Inputs:       b - b coefficients from the output of glmfit
%         covar_m - a covariate matrix where each column represents the
%                   data of a certain covariate at each timestep
%
% Outputs: lambda - conditional intensity at each timestep
%
% Function by: Lim Jing Xuan
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function lambda = gen_lambda(b, covar_m)

% preallocate for speed
loglambda = zeros(size(covar_m,1),1);

% add the fudge factor
loglambda = loglambda + b(1);

for i=1:size(covar_m,2)
   loglambda = loglambda + b(i+1)*covar_m(:,i);
end

lambda = exp(loglambda);
end