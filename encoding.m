%% start anew
clearvars; % clear previous variables
close all; % close previous plots

%% load data & generate variables
load('train.mat')
[ISI,XLocAtSpikes,YLocAtSpikes] =...
    gen_spike_stat(xN,yN,spikes_binned);

% wait bar prep
formatOut = 'yymmdd-HHMMSS';
date = datestr(now,formatOut);

%% plotting raw data
% spike times vs location
plot_spiking_positions(xN,yN,XLocAtSpikes,YLocAtSpikes,'subplot');
saveas(gcf, 'spike_pos.png')

% ISIs
ISI_threshold = 600;
ISIs = plot_ISIs(spikes_binned,ISI_threshold,2);
saveas(gcf, 'isi.png')

%% generate new covariates
[Vx,Vy,phi,r] = generate_new_variables(xN,yN,1000);  % raw data

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

% parameters
neurons = 1:10;
num_model = 3;
[x_new,y_new] = meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

% prepare storage
spike = cell(length(neurons),num_model);
covar = cell(length(neurons),num_model);
covar_grid = cell(length(neurons),num_model);
b = cell(length(neurons),num_model);
dev = cell(length(neurons),num_model);
stats = cell(length(neurons),num_model);
lambda = cell(length(neurons),num_model);
lambda_grid = cell(length(neurons),num_model);

%% iterate through neurons of interest

for i = neurons
    disp(['Working on neuron ' num2str(i) ' ...'])
    
    %% variables
    spikes = spikes_binned(:,i);  % spikes of the relevant neuron
    

    %% DEFINE MODELS
    % Model 1: neuron 6
    hist = [4:15 96:109 140:146];
    [spike{i,1},covar{i,1}] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,vyN,phi);
    [b{i,1},dev{i,1},stats{i,1}] = glmfit(covar{i,1},spike{i,1},'poisson');
    lambda{i,1} = gen_lambda(b{i,1},covar{i,1});
    % plot lambda as a function of X and Y position
    covar_grid{i,1} = [x_new,y_new,x_new.^2,y_new.^2,x_new.*y_new]';
    lambda_grid{i,1} = exp(b{i,1}(1) + ...
                           b{i,1}(2)*x_new + ...
                           b{i,1}(3)*y_new + ...
                           b{i,1}(4)*x_new.^2 + ...
                           b{i,1}(5)*y_new.^2 + ...
                           b{i,1}(6)*x_new.*y_new);
    lambda_grid{i,1}(find(x_new.^2 + y_new.^2 > 1)) = nan;
    
    % Model 2: unimodal place cells 1-5
    hist = [3:29 88:138];
    [spike{i,2},covar{i,2}] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,vxN,r,phi.^2);
    [b{i,2},dev{i,2},stats{i,2}] = glmfit(covar{i,2},spike{i,2},'poisson');
    lambda{i,2} = gen_lambda(b{i,2},covar{i,2});
    % plot lambda as a function of X and Y position
    covar_grid{i,2} = [x_new,y_new,x_new.^2,y_new.^2,x_new.*y_new]';
    lambda_grid{i,2} = exp(b{i,2}(1) + ...
                           b{i,2}(2)*x_new + ...
                           b{i,2}(3)*y_new + ...
                           b{i,2}(4)*x_new.^2 + ...
                           b{i,2}(5)*y_new.^2 + ...
                           b{i,2}(6)*x_new.*y_new);
    lambda_grid{i,2}(find(x_new.^2 + y_new.^2 > 1)) = nan;
    
    % Model 3: multimodal place cell 7-10
    hist = [4:30 96:146];
    [spike{i,3},covar{i,3}] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN,r,phi);
    [b{i,3},dev{i,3},stats{i,3}] = glmfit(covar{i,3},spike{i,3},'poisson');
    lambda{i,3} = gen_lambda(b{i,3},covar{i,3});
    % plot lambda as a function of X and Y position
    covar_grid{i,3} = [x_new,y_new,x_new.^2,y_new.^2,x_new.*y_new]';
    lambda_grid{i,3} = exp(b{i,3}(1) + ...
                           b{i,3}(2)*x_new + ...
                           b{i,3}(3)*y_new + ...
                           b{i,3}(4)*x_new.^2 + ...
                           b{i,3}(5)*y_new.^2 + ...
                           b{i,3}(6)*x_new.*y_new);
    lambda_grid{i,3}(find(x_new.^2 + y_new.^2 > 1)) = nan;
    
    waitbar(i/length(neurons),h);
    
end

for i = neurons
    disp(['Plotting neuron ' num2str(i) ' ...'])
    
    %% EVALUATE MODELS
    
    figure('units','normalized','outerposition',[0 0.035 1 0.92]);
    suptitle(['Cell ' num2str(i)]);

    % Model subplots w/beta subplots below
    for j = 1:num_model
        
        % plot lambda as a function of position
        subplot(2,num_model,j); hold on;
        %% TODO
        % the lambdas could be the same if they correspond to the same
        % variate (check)
        
        h_mesh = mesh(x_new,y_new,lambda_grid{i,j},'AlphaData',0);
        get(h_mesh,'AlphaData');
        set(h_mesh,'AlphaData',0);
        hold on;
        plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
        xlabel('x position [m]'); ylabel('y position [m]');
        
        % plot beta
        subplot(2,num_model,j+3); hold on;
        errorbar(b{i,j},2*stats{i,j}.se);
        xticks(1:length(b{i,j}));
        xlim([0 length(b{i,j})+1]);
        xticks(0:10:length(b{i,j}));
        xlabel('\beta number'); ylabel('\beta value');
        
        set(gca,'FontSize',16)
        saveas(gcf, [date '-beta_' num2str(i) '.png'])
        
    end
    
    %% plot KS plots for all three models
    ks_spikes{1} = spike{i,1};
    ks_spikes{2} = spike{i,2};
    ks_spikes{3} = spike{i,3};
    ks_lambda{1} = lambda{i,1};
    ks_lambda{2} = lambda{i,2};
    ks_lambda{3} = lambda{i,3};
    ks_dev{1} = dev{i,1};
    ks_dev{2} = dev{i,2};
    ks_dev{3} = dev{i,3};
    ks_b{1} = b{i,1};
    ks_b{2} = b{i,2};
    ks_b{3} = b{i,3};
    plot_ks(ks_spikes,ks_lambda,ks_b,ks_dev);
    cur_title = get(gca, 'Title');
    title([cur_title.String ': neuron ' num2str(i)]);
    
    set(gca,'FontSize',16)
    saveas(gcf, [date '-KS-neuron_' num2str(i) '.png'])
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
