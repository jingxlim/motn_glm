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

%% downsample data
ds_rate = 50;  % Hz
[xN_ds,yN_ds,spikes_binned_ds] = downsample(xN,yN,spikes_binned,ds_rate);

%% generate new covariates
[Vx,Vy,dir,r] = generate_new_variables(xN,yN,spikes_binned,1000);  % raw data
[Vx_ds,Vy_ds,dir_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,spikes_binned_ds,ds_rate);

clear spikess
clear lambdaEst

h = waitbar(0,'Please wait...');
neurons = [3 5 6 7 8 9 10];

for i=neurons  % iterate through all the neurons
    
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
    
%     spikes_ds = spikes_binned_ds(:,i);
%     % Model 3: linear + quadratic (downsampled)
%     covar_m3 = [xN_ds yN_ds xN_ds.^2 yN_ds.^2];
%     [b3,dev3,stats3] = glmfit(covar_m3,spikes_ds,'poisson');
%     lambdaEst{3} = gen_lambda(b3,covar_m3);
%     spikess{3} = spikes_ds;
%     
%     % Model 4: linear + quadratic + integrate (downsampled)
%     covar_m4 = [xN_ds yN_ds xN_ds.^2 yN_ds.^2 xN_ds.*yN_ds];
%     [b4,dev4,stats4] = glmfit(covar_m4,spikes_ds,'poisson');
%     lambdaEst{4} = gen_lambda(b4,covar_m4);
%     spikess{4} = spikes_ds;    
    
    % Model 3: linear + quadratic + integrate + history dependence
    hist = 3:180;
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
    lambda3 = exp(b3(1) + b3(2)*x_new + b3(3)*y_new + b3(4)*x_new.^2 +...
                  b3(5)*y_new.^2 + b3(6)*x_new.*y_new);
    lambda3(find(x_new.^2 + y_new.^2 > 1))=nan;  % what is this line doing?

    % plot lambda as a function position over this grid
    h_mesh = mesh(x_new,y_new,lambda3,'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
    xlabel('x position [m]'); ylabel('y position [m]');
    
    % plot betas
    subplot(1,2,2); hold on;
%     for n=1:numel(b3)
%         plot(n,b3(n),'*','DisplayName',num2str(stats3.p(n)));
%     end
%     legend('show','Location','bestoutside')
    errorbar(b3,2*stats3.se);
    xticks(1:length(b3));
    xlim([0 length(b3)+1]);
    xlabel('\beta number'); ylabel('\beta value');
    saveas(gcf, [date '-betas-neuron_' num2str(i) '.png'])
    save([date '-glm_out-neuron_' num2str(i) '.mat'],'b3','dev3','stats3')
    
    % plot KS plots for all three models
    plot_ks(spikess,lambdaEst);
    cur_title = get(gca, 'Title');
    title([cur_title.String ': neuron ' num2str(i)]);
    saveas(gcf, [date '-KS-neuron_' num2str(i) '.png'])
    
    waitbar(i/size(neurons,2),h);
end