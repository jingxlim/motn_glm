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

for i=1:size(covar_m1,2)
   loglambda = b(i+1).*covar_m(i);
end

log_lambda = log_lambda + b(0);

lambda = exp(log_lambda);
end