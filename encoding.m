%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models of the Neuron
% Project 2
% Group: Jingxuan Lim, Simon Orozco, Seony Han
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start anew
clearvars; % clear previous variables
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
clear spikess
clear lambdaEst

h = waitbar(0,'Please wait...');
for i=1:size(spikes_binned,2)  % iterate through all the neurons
    
    disp(['Working on neuron ' num2str(i) ' ...'])
    spikes = spikes_binned(:,i);
    
    % DEFINE MODELS
    % Model 1: linear + quadratic
    covar_m1 = [xN yN xN.^2 yN.^2];
    [b1,dev1,stats1] = glmfit(covar_m1,spikes,'poisson');
    lambdaEst{1} = gen_lambda(b1,covar_m1);
    spikess{1} = spikes;
    
    % Model 2: linear + quadratic + integrate
    covar_m2 = [xN yN xN.^2 yN.^2 xN.*yN];
    [b2,dev2,stats2] = glmfit(covar_m2,spikes,'poisson');
    lambdaEst{2} = gen_lambda(b2,covar_m2);
    spikess{2} = spikes;
    
    % Model 3: linear + quadratic + integrate + history dependence
    hist = 1:120;
    [spikes_m3,covar_m3] = hist_dep(hist,spikes,xN,yN,xN.^2,yN.^2,xN.*yN);
    [b3,dev3,stats3] = glmfit(covar_m3,spikes_m3,'poisson');
    lambdaEst{3} = gen_lambda(b3,covar_m3);
    spikess{3} = spikes_m3;
    
    % EVALUATE MODELS
    figure(); clf; hold on;
    set(gcf,'units','points','position',[100,100,1000,400])

    subplot(1,2,1); hold on;
    [x_new,y_new] = meshgrid(-1:.1:1);
    y_new = flipud(y_new);
    x_new = fliplr(x_new);

    % compute lambda for each point on this grid using the GLM model
    lambda3 = lambdaEst{3};
    lambda3(x_new.^2 + y_new.^2 > 1)=nan;

    % plot lambda as a function position over this grid
    h_mesh = mesh(x_new,y_new,lambda3,'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
    xlabel('x position [m]'); ylabel('y position [m]');
    
    % plot betas
    subplot(1,2,2); hold on;
    for n=1:numel(b3)
        plot(n,b3(n),'*','DisplayName',num2str(stats3.p(n)));
    end
    legend('show','Location','bestoutside')
    errorbar(b3,2*stats3.se);
    xticks(1:length(b3));
    xlim([0 length(b3)+1]);
    xlabel('\beta number'); ylabel('\beta value');
    saveas(gcf, ['betas-neuron_' num2str(i) '.png'])
    
    % plot KS plots for all three models
    plot_ks(spikess,lambdaEst);
    cur_title = get(gca, 'Title');
    title([cur_title.String ': neuron ' num2str(i)]);
    saveas(gcf, ['KS-neuron_' num2str(i) '.png'])
    
    waitbar(i/size(spikes_binned,2),h)
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

%% validate covariate choices

