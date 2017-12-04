%% start anew
clearvars; % clear previous variables
close all; % close previous plots

%% load data & generate variables

load('train.mat')
[ISI,XLocAtSpikes,YLocAtSpikes] =...
    gen_spike_stat(xN,yN,spikes_binned);


%%%% notes %%%%
% more variables could be useful but also means more to present

formatOut = 'yymmdd';
date = datestr(now,formatOut);

%% plotting raw data
% plot spike times as a function of location
plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes,'subplot');

% plot ISIs
ISI_threshold = 600;
ISIs = plot_ISIs(spikes_binned,ISI_threshold,2);

% plot_spiking_velocities maybe?
% plot_spiking_directions maybe?

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
% encode_place
% encode_grid
% encode_mm
% ^ idk if those are functions. read below
clear spikess
clear lambdaEst

h = waitbar(0,'Please wait...');

neurons = 10; models = 3;

for i = 1:neurons  % iterate through all the neurons
    % variables
    b = cell(1,models);
    stats = cell(1,models);
    
    disp(['Working on neuron ' num2str(i) ' ...'])
    spikes = spikes_binned(:,i);
    
    % DEFINE MODELS
    % Model 1: neuron 6
    hist = [7:18 99:112 143:149];
    [spikes_m1,covar_m1] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,phi);
    [b{1},dev1,stats{1}] = glmfit(covar_m1,spikes_m1,'poisson');
    lambdaEst{1} = gen_lambda(b{1},covar_m1);
    spikess{1} = spikes_m1;
    
    % Model 2: unimodal place cells 1-5
    hist = [3:30 91:141];
    [spikes_m2,covar_m2] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,phi.^2);
    [b{2},dev2,stats{2}] = glmfit(covar_m2,spikes_m2,'poisson');
    lambdaEst{2} = gen_lambda(b{2},covar_m2);
    spikess{2} = spikes_m2;
    
    % Model 3: multimodal place cell 7-10
    hist = [7:33 99:149];
    [spikes_m3,covar_m3] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,phi);
    [b{3},dev3,stats{3}] = glmfit(covar_m3,spikes_m3,'poisson');
    lambdaEst{3} = gen_lambda(b{3},covar_m3);
    spikess{3} = spikes_m3;
    
    % EVALUATE MODELS
    figure('Name',['Cell ' num2str(i)]);
    [x_new,y_new] = meshgrid(-1:.1:1);

    % Model subplots w/beta subplots below
    for j = 1:models
        % model
        subplot(models,2,j); hold on;
        h_mesh = mesh(x_new,y_new,lambdaEst{j},'AlphaData',0);
        plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
        xlabel('x position [m]'); ylabel('y position [m]');
        % beta
        subplot(models,2,j+models); hold on;
        errorbar(b{j},2*stats{j}.se);
        xticks(1:length(b{j}));
        xlim([0 length(b{j})+1]);
        xlabel('\beta number'); ylabel('\beta value');
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
