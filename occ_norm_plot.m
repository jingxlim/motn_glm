%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% occ_norm_plot(spikes_binned,vxN,vyN,phi,r)
% 
% Outputs: none (produces figures)
% Inputs:
%   spikes_binned - spikes for all neurons (change to individual later)
%             vxN - x-direction velocities of rat
%             vyN - y-direction velocities of rat
%             phi - movement directions of rat
%               r - movement speed of rat
% 
% This function creates occupancy normalized plots for all neurons using
% the velocities, movement directions, and speed of the animal. Different
% models are tested on each parameter vxN, vyN, phi, & r. The R-S plots of
% each model is plotted on a single figure for each parameter.
% 
% When using this function from Command Prompt:
%   load('train.mat');
%   [vxN,vyN,phi,r] = generate_new_variables(xN,yN,1000);
% 
% Adapted from glm_part2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function occ_norm_plot(spikes_binned,vxN,vyN,phi,r)
% models tested
m = 2; % number of models tested
p = 4; % number of parameters involved

% more variable stuff for efficiency-ish
params{1} = vxN;
params{2} = vyN;
params{3} = phi;
params{4} = r;
names{1} = 'vxN';
names{2} = 'vyN';
names{3} = 'phi';
names{4} = 'r';
m_names{1} = 'linear';
m_names{2} = 'quad';
figs = cell(1,m*p);


% testing models
% each parameter
for i = 1:p
    % add/remove lines here based on # of models tested
    figs{m*i-1} = figure('Name',['ONP: ' names{i} ' unimodal']);
    figs{m*i}   = figure('Name',['ONP: ' names{i} ' multimodal']);
    % each neuron
    for j = 1:10
        test_mods(spikes_binned,params{i},names{i},m_names,j,m,figs,i*m-1);
    end
    saveas(gcf, ['ONP_' num2str(names{ceil(i/m)}) '_' num2str(m_names{i-2*floor((i-1)/m)}) '.png'])
end

%% saving figures
% for i = 1:8
%     figure(figs{i})
%     saveas(gcf, ['ONP_' num2str(names{ceil(i/m)}) '_' num2str(m_names{i-2*floor((i-1)/m)}) '.png'])
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test(spikes_binned,data,name,n)
% 
% Outputs: none (figures only)
% Inputs:
%   spikes_binned - spikes for all neurons (change to individual later)
%   data - parameter used
%   name - String name of data used
%   n - neuron #
% 
% This function tests 6 different covariate models using the parameter data
% given on the n neuron and calculates error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_mods(spikes_binned,data,name,m_names,n,m,figs,fig_id)
% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:double(max(spikes_binned))
    ind = [ind; find(spikes_binned(:,n) >= numspikes)];
end;

% variables
b = cell(1,m);
cov = cell(1,m);
data_x = min(data):range(data)/80:max(data);

%% histogram
for i = 0:m-1
    figure(figs{fig_id+i})
    subplot(2,5,n)
    bar(data_x,hist(data(ind),data_x)./hist(data,data_x));
    title(['Neuron ' num2str(n)]) % ', ' name ', ' m_names{i+1}])
    if n==8
        xlabel(name);
    end
    if mod(n,5)==1
        ylabel('normalized spike counts');
    end
end

%% models
% model 1: linear
cov{1} = [data];
b{1} = glmfit(cov{1},spikes_binned(:,n),'poisson');
% figure
figure(figs{fig_id})
hold on;
subplot(2,5,n)
plot(data_x,exp(b{1}(1)+b{1}(2)*data_x),'r','LineWidth',1.5);
set(gca,'fontsize',20)

% model 2: quad
cov{2} = [data.^2];
b{2} = glmfit(cov{2},spikes_binned(:,n),'poisson');
% figure
figure(figs{fig_id+1})
hold on;
subplot(2,5,n)
plot(data_x,exp(b{2}(1)+b{2}(2)*data_x.^2),'r','LineWidth',1.5);
set(gca,'fontsize',20)

%% error calculation: creates m*10 figures
% find_ks(b,cov,spikes_binned,name,n,m);

%% other models
% % model 3
% cov{3} = [data data.^2];
% b{3} = glmfit(cov{3},spikes_binned(:,n),'poisson');
% hold on;
% plot(data_x,exp(b{3}(1)+b{3}(2)*data_x+b{3}(3)*data_x.^2),'c');
% 
% % model 4
% cov{4} = [data.^3];
% b{4} = glmfit(cov{4},spikes_binned(:,n),'poisson');
% hold on;
% plot(data_x,exp(b{4}(1)+b{4}(2)*data_x.^3),'m');
% 
% % model 5
% cov{5} = [data data.^3];
% b{5} = glmfit(cov{5},spikes_binned(:,n),'poisson');
% hold on;
% plot(data_x,exp(b{5}(1)+b{5}(2)*data_x+b{5}(3)*data_x.^3),'b');
% 
% % model 6
% cov{6} = [data data.^2 data.^3];
% b{6} = glmfit(cov{6},spikes_binned(:,n),'poisson');
% hold on;
% plot(data_x,exp(b{6}(1)+b{6}(2)*data_x+b{6}(3)*data_x.^2+b{6}(4)*data_x.^3),'y');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_ks(b,cov,spikes_binned,name,n,m)
% 
% Outputs: none (figures)
% Inputs:
%   b - vector of beta values
%   cov - vector of covariates
%   name - String name of data used
%   n - neuron #
%   m - # of models
% 
% This function creates the lambdaEst and spikes matrix necessary to call
% the plot_ks function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function find_ks(b,cov,spikes_binned,name,n,m)
% variables
lambdaEst = cell(1,m);
spikes = cell(1,m);

% loop to fill variables
for i = 1:m
    lambdaEst{i} = gen_lambda(b{i},cov{i});
    spikes{i} = spikes_binned;
end

% actually calculations
plot_ks(spikes,lambdaEst);
cur_title = get(gca, 'Title');
title([cur_title.String ': neuron ' num2str(n) ' ' name]);
% saveas(gcf, ['KS-neuron_' num2str(n) '.png'])

end