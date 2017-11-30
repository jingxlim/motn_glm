%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from glm_part2.m
% 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% delete when using as function
% load the rat trajectory and spiking data
clearvars
load('train.mat');
% downsample
ds_rate = 50;  % Hz
[xN_ds,yN_ds,spikes_binned_ds] = downsample(xN,yN,spikes_binned,ds_rate);
% variable generation
[Vx,Vy,phi,r] = generate_new_variables(xN,yN,spikes_binned,1000);  % raw data
[Vx_ds,Vy_ds,phi_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,spikes_binned_ds,ds_rate);
% this will be inputs to the function
n = 1; % neuron #

%% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:max(spikes_binned)
    ind = [ind; find(spikes_binned(:,n) >= numspikes)];
end;

figure('Name','Occupancy Normalized Histogram');
warning('off');

%% x-velocity
% Histogram of spiking to x-velocity
subplot(2,2,1);
    velocities = min(Vx):range(Vx)/80:max(Vx);
    bar(velocities,hist(Vx(ind),velocities)./hist(Vx,velocities));
    xlabel('x velocity');
    ylabel('normalized spike counts');

    % GLM
    b = glmfit(Vx,spikes_binned(1:end-1,n),'poisson');
    hold on;
    plot(velocities,exp(b(1)+b(2)*-hist(Vx(ind),velocities)./hist(Vx,velocities)),'r');
    
%% y-velocity
% Histogram of spiking to y-velocity 
subplot(2,2,2);
    velocities = min(Vy):range(Vy)/80:max(Vy);
    bar(velocities,hist(Vy(ind),velocities)./hist(Vy,velocities));
    xlabel('y velocity');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(Vy,spikes_binned(1:end-1,n),'poisson');
    hold on;
   plot(velocities,exp(b(1)+b(2)*hist(Vy(ind),velocities)./hist(Vy,velocities)),'r');

%% movement speed
% Histogram of spiking to movement speed 
subplot(2,2,3);
    rs = 0:range(r)/40:max(r);
    bar(rs,hist(r(ind),rs)./hist(r,rs));
    xlabel('speed');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(r,spikes_binned(1:end-1,n),'poisson');    
    hold on;
    plot(rs,exp(b(1)+b(2)*hist(r(ind),rs)./hist(r,rs)),'r');

%% movement direction
% Histogram of spiking to movement direction
subplot(2,2,4);
    phis = -pi:.1:pi;
    bar(phis,hist(phi(ind),phis)./hist(phi,phis));
    xlabel('direction');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(phi,spikes_binned(1:end-1,n),'poisson');    
    hold on;
    plot(phis,exp(b(1)+b(2)*hist(phi(ind),phis)./hist(phi,phis)),'r');
