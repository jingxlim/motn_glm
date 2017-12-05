%% start anew
clearvars; % clear previous variables
close all; % close previous plots

%% load data & generate variables
load('train.mat')
[ISI,XLocAtSpikes,YLocAtSpikes] =...
    gen_spike_stat(xN,yN,spikes_binned);
neurons = size(spikes_binned,2);

% wait bar prep
formatOut = 'yymmdd';
date = datestr(now,formatOut);

%% plotting raw data
% spike times vs location
plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes,'subplot');

% ISIs
ISI_threshold = 600;
ISIs = plot_ISIs(spikes_binned,ISI_threshold,2);

%% downsample data
ds_rate = 50;  % Hz
[xN_ds,yN_ds,spikes_binned_ds] = downsample(xN,yN,spikes_binned,ds_rate);

%% generate new covariates
[Vx,Vy,phi,r] = generate_new_variables(xN,yN,1000);  % raw data
[Vx_ds,Vy_ds,dir_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,ds_rate);

%% classifying cells
% This needs to be automated, and it turns out that is a hard problem (at
% least I struggled with it -Simon)

% place
% grid
% multimodal (?) (I think this is where the history dependence could prove
% to be a helpful covariate -Simon)

%% encoding data
clear spikess
clear lambdaEst

h = waitbar(0,'Please wait...');
models = 3;

% iterate through all the neurons
for i = 1:neurons
    disp(['Working on neuron ' num2str(i) ' ...'])
    
    %% variables
    spikes = spikes_binned(:,i);
    spikess = cell(1,models);
    covar = cell(1,models);
    [x_new,y_new] = meshgrid(-1:.1:1);

    
    % w/o hist dep (for lambda plot)
    b = cell(1,models);
    lambda = cell(1,models);
    
    % w/hist dep
    b_hist = cell(1,models);
    stats = cell(1,models);
    
    %% DEFINE MODELS
    % Model 1: neuron 6
    covar{1} = [xN,yN,xN.^2,yN.^2,xN.*yN,phi];
    hist = [7:18 99:112 143:149];
    lambda{1} = ones(length(xN),1);
    
    [spikess{1},covar_m1] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,phi);
    [b_hist{1},dev1,stats{1}] = glmfit(covar_m1,spikess{1},'poisson');
    lambdaEst{1} = gen_lambda(b_hist{1},covar_m1);
    % prep for lambda plot
    [b{1},~,~] = glmfit(covar{1},spikes,'poisson');
    for j = 1:length(covar{1})
        lambda{1} = b{1}(j)*covar{1}(:,j) + lambda{1};
    end
    lambda{1} = exp(lambda{1});
%     lambda{1} = exp(b{1}(1) + b{1}(2)*x_new + b{1}(3)*y_new + ...
%         b{1}(4)*x_new.^2 + b{1}(5)*y_new.^2 + b{1}(6)*x_new.*y_new + ...
%         b{1}(7)*phi);
    lambda{1}(find(x_new.^2 + y_new.^2 > 1)) = nan; % Simon what does this mean?
    
    % Model 2: unimodal place cells 1-5
    covar{2} = [ones(length(xN),1),xN,yN,xN.^2,yN.^2,xN.*yN,phi.^2];
    hist = [3:30 91:141];
    lambda{2} = [];
    
    [spikess{2},covar_m2] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,phi.^2);
    [b_hist{2},dev2,stats{2}] = glmfit(covar_m2,spikess{2},'poisson');
    lambdaEst{2} = gen_lambda(b_hist{2},covar_m2);
    % prep for lambda plot
    covar{2} = [ones(length(xN),1),xN,yN,xN.^2,yN.^2,xN.*yN,phi.^2];
    [b{2},~,~] = glmfit(covar{2},spikes,'poisson');
    for j = 1:length(covar{2})
        lambda{2} = b{2}(j)*covar{2}(:,j) + lambda{2};
    end
    lambda{2} = exp(lambda{2});
    lambda{2}(find(x_new.^2 + y_new.^2 > 1)) = nan; % Simon what does this mean?

    
    % Model 3: multimodal place cell 7-10
    covar{3} = [ones(length(xN),1),xN,yN,xN.^2,yN.^2,xN.*yN,phi];
    hist = [7:33 99:149];
    lambda{3} = [];

    [spikess{3},covar_m3] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,phi);
    [b_hist{3},dev3,stats{3}] = glmfit(covar_m3,spikess{3},'poisson');
    lambdaEst{3} = gen_lambda(b_hist{3},covar_m3);
    % prep for lambda plot
    [b{3},~,~] = glmfit(covar{3},spikes,'poisson');
    for j = 1:length(covar{2})
        lambda{3} = b{3}(j)*covar{3}(:,j) + lambda{3};
    end
    lambda{3} = exp(lambda{3});
    lambda{3}(find(x_new.^2 + y_new.^2 > 1)) = nan; % Simon what does this mean?

    %% EVALUATE MODELS
    figure('Name',['Cell ' num2str(i)]);

    % Model subplots w/beta subplots below
    for j = 1:models
        % model
        subplot(models,2,j); hold on;
        %% TODO
        % create lambda w/o history for each model
        % the lambdas could be the same if they correspond to the same
        % variate (check)
        % use x_new and y_new to calculate the lambda --> same size
        
        %%
        h_mesh = mesh(x_new,y_new,lambdaEst{j},'AlphaData',0);
        plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
        xlabel('x position [m]'); ylabel('y position [m]');
        % beta
        subplot(models,2,j+models); hold on;
        errorbar(b_hist{j},2*stats{j}.se);
        xticks(1:length(b_hist{j}));
        xlim([0 length(b_hist{j})+1]);
        xlabel('\beta number'); ylabel('\beta value');
        % k-s
        
        % 
    end
    
    %% was commented out before
%     for n=1:numel(b3)
%         plot(n,b3(n),'*','DisplayName',num2str(stats3.p(n)));
%     end
%     legend('show','Location','bestoutside')

    %% put into loop above
%     errorbar(b3,2*stats3.se);
%     xticks(1:length(b3));
%     xlim([0 length(b3)+1]);
%     xlabel('\beta number'); ylabel('\beta value');
%     saveas(gcf, [date '-betas-neuron_' num2str(i) '.png'])
%     save([date '-glm_out-neuron_' num2str(i) '.mat'],'b3','dev3','stats3')
    
    %% plot KS plots for all three models
    plot_ks(spikess,lambdaEst);
    cur_title = get(gca, 'Title');
    title([cur_title.String ': neuron ' num2str(i)]);
    saveas(gcf, [date '-KS-neuron_' num2str(i) '.png'])
    
    waitbar(i/length(neurons),h);
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
