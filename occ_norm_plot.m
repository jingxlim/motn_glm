%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from glm_part2.m
% 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% delete when using as function
% load the rat trajectory and spiking data
% clearvars
% load('train.mat');
% downsample
ds_rate = 50;  % Hz
% [xN_ds,yN_ds,spikes_binned_ds] = downsample(xN,yN,spikes_binned,ds_rate);
% variable generation
% [vxN,vyN,phi,r] = generate_new_variables(xN,yN,spikes_binned,1000);  % raw data
% [Vx_ds,Vy_ds,phi_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,spikes_binned_ds,ds_rate);
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
%     velocities = min(vxN):range(vxN)/80:max(vxN);
    velocities = -.04:.001:.04;
    bar(velocities,hist(vxN(ind),velocities)./hist(vxN,velocities));
    xlabel('x velocity');
    ylabel('normalized spike counts');

    % GLM
    b = glmfit(vxN,spikes_binned(:,n),'poisson');
    hold on;
%     plot(velocities,exp(b(1)+b(2)*hist(vxN(ind),velocities)./hist(vxN,velocities)),'r');
%     plot(velocities,exp(b(1)+b(2)*vxN(ind)),'r');
    plot(velocities,exp(b(1)+b(2)*velocities),'r');

%% y-velocity
% Histogram of spiking to y-velocity 
subplot(2,2,2);
%     velocities = min(vyN):range(vyN)/80:max(vyN);
    bar(velocities,hist(vyN(ind),velocities)./hist(vyN,velocities));
    xlabel('y velocity');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(vyN,spikes_binned(:,n),'poisson');
    hold on;
%     plot(velocities,exp(b(1)+b(2)*hist(vyN(ind),velocities)./hist(vyN,velocities)),'r');
    plot(velocities,exp(b(1)+b(2)*velocities),'r');

%% movement speed
% Histogram of spiking to movement speed 
subplot(2,2,3);
%     rs = 0:range(r)/40:max(r);
    rs = 0:.001:.04;
    bar(rs,hist(r(ind),rs)./hist(r,rs));
    xlabel('speed');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(r,spikes_binned(:,n),'poisson');    
    hold on;
%     plot(rs,exp(b(1)+b(2)*hist(r(ind),rs)./hist(r,rs)),'r');
    plot(rs,exp(b(1)+b(2)*rs),'r');

%% movement direction
% Histogram of spiking to movement direction
subplot(2,2,4);
    phis = -pi:.1:pi;
    bar(phis,hist(phi(ind),phis)./hist(phi,phis));
    xlabel('direction');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(phi,spikes_binned(:,n),'poisson');    
    hold on;
%     plot(phis,exp(b(1)+b(2)*hist(phi(ind),phis)./hist(phi,phis)),'r');
    plot(phis,exp(b(1)+b(2)*phis),'r');
