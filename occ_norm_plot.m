%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% occ_norm_plot(spikes_binned,vxN,vyN,phi,r,n)
% 
% Outputs: none (produces figures)
% Inputs:
%   spikes_binned - spikes for all neurons (change to individual later)
%             vxN - x-direction velocities of rat
%             vyN - y-direction velocities of rat
%             phi - movement directions of rat
%               r - movement speed of rat
%               n - neuron #
% 
% This function creates occupancy normalized plots for neuron n using the
% velocities, movement directions, and speed of the animal. Different
% models are tested on each parameter vxN, vyN, phi, & r. The R-S plots of
% each model is plotted on a single figure for each parameter.
% 
% When using this function from Command Prompt:
%   load('train.mat');
%   [vxN,vyN,phi,r] = generate_new_variables(xN,yN,spikes_binned,1000);
% 
% Notes to self
%   Need to try history dependence
%   do downsampling later
%   ds_rate as input (Hz)
% 
% Adapted from glm_part2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function occ_norm_plot(spikes_binned,vxN,vyN,phi,r,n)%,ds_rate)
%% downsampling stuff for later
% [xN_ds,yN_ds,spikes_binned_ds] = downsample(xN,yN,spikes_binned,ds_rate);
% [Vx_ds,Vy_ds,phi_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,spikes_binned_ds,ds_rate);

%% modify variables to work in functions
vxN = [vxN;0];
vyN = [vyN;0];
r = [r;0];
phi = [phi;0];

%% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:max(spikes_binned)
    ind = [ind; find(spikes_binned(:,n) >= numspikes)];
end;

% no history dependence
test(spikes_binned,vxN,'vxN',ind,n);
test(spikes_binned,vyN,'yxN',ind,n);
test(spikes_binned,r,  'r',  ind,n);
test(spikes_binned,phi,'phi',ind,n);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%   data - parameter used
%   name - String name of data used
%   n - neuron #
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test(spikes_binned,data,name,ind,n)
figure('Name','Occupancy Normalized Plot')
data_x = min(data):range(data)/80:max(data);
bar(data_x,hist(data(ind),data_x)./hist(data,data_x));
xlabel(name);
ylabel('normalized spike counts');

% model 1
cov1 = [data];
b1 = glmfit(cov1,spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b1(1)+b1(2)*data_x),'r');

% model 2
cov2 = [data.^2];
b2 = glmfit(cov2,spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b2(1)+b2(2)*data_x.^2),'g');

% model 3
cov3 = [data data.^2];
b3 = glmfit(cov3,spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b3(1)+b3(2)*data_x+b3(3)*data_x.^2),'c');

% model 4
cov4 = [data.^3];
b4 = glmfit(cov4,spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b4(1)+b4(2)*data_x.^3),'m');

% model 5
cov5 = [data data.^3];
b5 = glmfit(cov5,spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b5(1)+b5(2)*data_x+b5(3)*data_x.^3),'b');

% model 6
cov6 = [data data.^2 data.^3];
b6 = glmfit(cov6,spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b6(1)+b6(2)*data_x+b6(3)*data_x.^2+b6(4)*data_x.^3),'y');

% legend
title(['ONP for Neuron ' num2str(n) ' ' name])
legend('velocities','lin','quad','lin+quad','cub','li+cub','lin+quad+cub')
set(gca, 'fontsize',14)

% error calculation -> make this into loop later w/b&cov as cell arrays
lambdaEst{1} = gen_lambda(b1,cov1);
lambdaEst{2} = gen_lambda(b2,cov2);
lambdaEst{3} = gen_lambda(b3,cov3);
lambdaEst{4} = gen_lambda(b4,cov4);
lambdaEst{5} = gen_lambda(b5,cov5);
lambdaEst{6} = gen_lambda(b6,cov6);

spikes = cell(1,numel(lambdaEst));
for i = 1:numel(lambdaEst)
    spikes{i} = spikes_binned;
end

plot_ks(spikes,lambdaEst);
cur_title = get(gca, 'Title');
title([cur_title.String ': neuron ' num2str(n) ' ' name]);
% saveas(gcf, ['KS-neuron_' num2str(n) '.png'])
end