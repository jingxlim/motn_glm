
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decoder.m
% [x_predicted, y_predicted] = decoder(spikes_binned);
%
% This function takes a matrix of spike trains from n neurons and uses our
% previously established model to predict X and Y location from it.
%
% Code by: Simon Orozco
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_predicted, y_predicted] = decoder(spikes_binned)

T = size(spikes_binned,1);
n = size(spikes_binned,2);
hist = [3:29 88:138];
num_cov = 8+length(hist);

train_data = load('train.mat');
Ntrain = size(train_data.spikes_binned,2);

[Vx,~,phi,r] = generate_new_variables(train_data.xN,train_data.yN,1000);  % raw data
betas = NaN(num_cov+1,Ntrain);
lambdas = zeros(length(r),Ntrain);
for i = 1:Ntrain
    % Run the training model to get the beta values
    [spike,covar] = hist_dep(hist,train_data.spikes_binned(:,i),train_data.xN,train_data.yN,...
        train_data.xN.^2,train_data.yN.^2,train_data.xN.*train_data.yN,Vx,r,phi.^2);
    betas(:,i) = glmfit(covar,spike,'poisson');
    lambdas(max(hist)+1:end,i) = gen_lambda(betas(:,i),covar);
end

% Average across fits
betas = nanmean(betas,2);
lambda = nanmean(lambdas,2);

% Preallocate
covs = NaN(num_cov+1,T,n);
for i = 1:n
    % Calculate covariates from betas and conditional intensity
    for j = 1:length(betas)
        covs(j,:,i) = (betas(j)).\log(lambda);
        if j==1
            lambda = lambda - betas(j);
        else
            lambda = lambda - [zeros(max(hist),1); betas(j)*covar(:,j-1)];
        end
    end
end

% Calculate means of x and y using the covariates based on them
x_predicted = nanmean([nanmean(covs(1,:,:),3); nanmean(sqrt(covs(3,:,:)),3)]);
y_predicted = nanmean([nanmean(covs(2,:,:),3); nanmean(sqrt(covs(4,:,:)),3)]);

end


