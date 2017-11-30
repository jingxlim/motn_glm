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
[Vx,Vy,dir,r] = generate_new_variables(xN,yN,spikes_binned,1000);  % raw data
[Vx_ds,Vy_ds,dir_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,spikes_binned_ds,ds_rate);

%% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:max(spikes_binned)
    ind = [ind; find(spikes_binned >= numspikes)];
end;

figure('Name','Occupancy Normalized Histogram');
warning('off');

%% x-velocity
% Histogram of spiking to x-velocity
subplot(2,2,1);
    velocities = min(Vx):floor(range(Vx))/80:max(Vx);
    bar(velocities,hist(Vx(ind),velocities)./hist(Vx,velocities));
    xlabel('x velocity');
    ylabel('normalized spike counts');

    % GLM
    b = glmfit(vxN,spikes_binned,'poisson');
    hold on;
    plot(velocities,exp(b(1)+b(2)*-hist(Vx(ind),velocities)./hist(Vx,velocities)),'r');
    
%% y-velocity
% Histogram of spiking to y-velocity 
subplot(2,2,2);
    bar(velocities,hist(vyN(ind),velocities)./hist(vyN,velocities));
    xlabel('y velocity');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(vyN,spikes_binned,'poisson');
    hold on;
   plot(velocities,exp(b(1)+b(2)*hist(vyN(ind),velocities)./hist(vyN,velocities)),'r');

%% movement speed
% Histogram of spiking to movement speed 
subplot(2,2,3);
    rs = 0:.001:.04;
    bar(rs,hist(r(ind),rs)./hist(r,rs));
    xlabel('speed');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(r.^(0.5),spikes_binned,'poisson');    
    hold on;
    plot(rs,exp(b(1)+b(2)*(hist(r(ind),rs)./hist(r,rs)).^(0.5)),'r');

%% movement direction
% Histogram of spiking to movement direction
subplot(2,2,4);
    phis = -pi:.1:pi;
    bar(phis,hist(phi(ind),phis)./hist(phi,phis));
    xlabel('direction');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(phi,spikes_binned,'poisson');    
    hold on;
    plot(phis,exp(b(1)+b(2)*hist(phi(ind),phis)./hist(phi,phis)),'r');
