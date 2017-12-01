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
% Adapted from glm_part2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function occ_norm_plot(spikes_binned,vxN,vyN,phi,r,n)
%% modify variables to work in functions
vxN = [vxN;0];
vyN = [vyN;0];
r = [r;0];
phi = [phi;0];

% no history dependence
test(spikes_binned,vxN,'vxN',n);
test(spikes_binned,vyN,'yxN',n);
test(spikes_binned,r,  'r',  n);
test(spikes_binned,phi,'phi',n);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test(spikes_binned,data,name,ind,n)
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
function test(spikes_binned,data,name,n)
% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:max(spikes_binned)
    ind = [ind; find(spikes_binned(:,n) >= numspikes)];
end;

% figure prep
figure('Name','Occupancy Normalized Plot')
xlabel(name);
ylabel('normalized spike counts');

% variables
m = 6;
b = cell(1,m);
cov = cell(1,m);
data_x = min(data):range(data)/80:max(data);
bar(data_x,hist(data(ind),data_x)./hist(data,data_x));

% model 1
cov{1} = [data];
b{1} = glmfit(cov{1},spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b{1}(1)+b{1}(2)*data_x),'r');

% model 2
cov{2} = [data.^2];
b{2} = glmfit(cov{2},spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b{2}(1)+b{2}(2)*data_x.^2),'g');

% model 3
cov{3} = [data data.^2];
b{3} = glmfit(cov{3},spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b{3}(1)+b{3}(2)*data_x+b{3}(3)*data_x.^2),'c');

% model 4
cov{4} = [data.^3];
b{4} = glmfit(cov{4},spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b{4}(1)+b{4}(2)*data_x.^3),'m');

% model 5
cov{5} = [data data.^3];
b{5} = glmfit(cov{5},spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b{5}(1)+b{5}(2)*data_x+b{5}(3)*data_x.^3),'b');

% model 6
cov{6} = [data data.^2 data.^3];
b{6} = glmfit(cov{6},spikes_binned(:,n),'poisson');
hold on;
plot(data_x,exp(b{6}(1)+b{6}(2)*data_x+b{6}(3)*data_x.^2+b{6}(4)*data_x.^3),'y');

% legend
title(['ONP for Neuron ' num2str(n) ' ' name])
legend('velocities','lin','quad','lin+quad','cub','li+cub','lin+quad+cub')
set(gca, 'fontsize',14)

%% error calculation
lambdaEst = cell(1,m);
spikes = cell(1,m);

for i = 1:m
    lambdaEst{i} = gen_lambda(b{i},cov{i});
    spikes{i} = spikes_binned;
end

plot_ks(spikes,lambdaEst);
cur_title = get(gca, 'Title');
title([cur_title.String ': neuron ' num2str(n) ' ' name]);
% saveas(gcf, ['KS-neuron_' num2str(n) '.png'])
end
