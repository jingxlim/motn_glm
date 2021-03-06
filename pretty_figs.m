%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pretty_figs(spikes_binned,xN,yN,vxN,vyN,phi,r,n)
% 
% Outputs: none (produces figures)
% Inputs:
%   spikes_binned - spikes for all neurons (change to individual later)
%              xN - x-positions of rat
%              yN - y-positions of rat
%             vxN - x-direction velocities of rat
%             vyN - y-direction velocities of rat
%             phi - movement directions of rat
%               r - movement speed of rat
% 
% This function creates K-S plots for all neurons to compare different
% covariate models invovling the positions, velocities, movement
% directions, and speed of the animal. The K-S plots involving the same
% covariates are plotted on the same figure, which each figure displaying
% either only unimodal (neurons 1-5) or multimodal (neurons 6-10) cells.
% 
% When using this function from Command Prompt:
%   load('train.mat');
%   [vxN,vyN,phi,r] = generate_new_variables(xN,yN,1000);
% 
% Adapted from glm_part2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pretty_figs(spikes_binned,xN,yN,vxN,vyN,phi,r)
%% variables
% "global" variables
% m = 8; % number of models tested
p = 4; % number of parameters involved
% figs = cell(p,m/p);

% base: always included in covariates
base = [xN yN xN.^2 yN.^2 xN.*yN];

% covs: what we're testing
covs{1,1} = vxN;    covs{1,2} = vxN.^2;
covs{2,1} = vyN;    covs{2,2} = vyN.^2;
covs{3,1} = phi;    covs{3,2} = phi.^2;
covs{4,1} = r;      covs{4,2} = r.^2;

% for labeling purposes
p_names{1} = 'vxN'; p_names{2} = 'vyN';
p_names{3} = 'phi'; p_names{4} = 'r';
p_names{5} = 'base';
% modals{1} = 'uni';  modals{2} = 'multi';

%% testing
% each parameter in params
for i = 1:p
    % add/remove lines here based on # of models tested
%     figs{i,1} = figure('Name',['ONP: ' p_names{i} ' unimodal']);
%     figs{i,2} = figure('Name',['ONP: ' p_names{i} ' multimodal']);
    temp_covs{1} = covs{i,1};
    temp_covs{2} = covs{i,2};
    % each neuron
    for j = 1:10
        test_models(spikes_binned(:,j),base,temp_covs,p_names{i},j);%,figs);
    end
end

%% saving figures
% for i = 1:8
%     figure(figs{i})
%     saveas(gcf, ['ONP_' num2str(names{ceil(i/m)}) '_' num2str(m_names{i-2*floor((i-1)/m)}) '.png'])
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_models(spikes,base,covs,modals,n,figs)
% 
% Outputs: none (figures only)
% Inputs:
%   spikes_binned - spikes for neuron n
%   base - the covariates common in all models tested
%   covs - the covariates we are testing
%   modals - {'uni' 'multi'}
%   n - neuron #
%   p - parameter # (1=vxN, 2=vyN, 3=phi, 4=r)
%   figs - cell array of all the figures
% 
% This function tests 6 different covariate models using the parameter data
% given on the n neuron and calculates error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_models(spikes,base,covs,p_name,n)
% variables
m = numel(covs)+1;  % number of models testing = lin+quad+base
m_add = 2;          % if testing additional models
m_tot = m + m_add;
b = cell(1,m_tot);
cov = cell(1,m_tot);

%% models
% base
cov{1} = base;
b{1} = glmfit(cov{1},spikes,'poisson');

% additional covariates
for i = 2:m
    % glm
    cov{i} = [base covs{i-1}];
    b{i} = glmfit(cov{i},spikes,'poisson');    
end

% additional models (variations of base)
for i = 1:m_add
    ind = m+i;
    cov{ind} = base(:,i+1:end);
    b{ind} = glmfit(cov{ind},spikes,'poisson');
end
%% error calculation: creates m*10 figures
find_ks3(b,cov,spikes,p_name,n,m_tot);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_ks(b,cov,spikes_binned,name,n,m)
% 
% Outputs: none (figures)
% Inputs:
%   b - vector of beta values
%   cov - vector of covariates
%   spikes - spikes of neuron n
%   figs - ALL THE FIGURES MWAHAHAHAHAA (all the neurons & parameters)
%   mode - # representing modality (1=uni, 2=multi)
%   n - neuron #
%   p - parameter #
%   m - # of models
% 
% This function creates the lambdaEst and spikes matrix necessary to call
% the plot_ks function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function find_ks3(b,cov,spikes,p_name,n,m)
% variables
lambdaEst = cell(1,m);
spikes_all = cell(1,m);

% each model: base, linear quad
for i = 1:m
    lambdaEst{i} = gen_lambda(b{i},cov{i});
    spikes_all{i} = spikes;
end

% prep the plot
% figure(figs{p,p_name})
% subplot(1,5,n-5*(p_name-1))
% hold on;

% actually calculations
name = ['KS: Neuron ' num2str(n) ' w/' p_name];
figure('Name', name)
anyplot_ks(spikes_all,lambdaEst);
title(name);
set(gca,'fontsize',20)

% save figures
saveas(gcf, ['KS_Neuron_' num2str(n) '_' p_name '.png'])


end