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
% [xN_ds,yN_ds,spikes_binned_ds] = downsample(xN,yN,spikes_binned,ds_rate);
% variable generation
[vxN,vyN,phi,r] = generate_new_variables(xN,yN,spikes_binned,1000);  % raw data
% [Vx_ds,Vy_ds,phi_ds,r_ds] = generate_new_variables(xN_ds,yN_ds,spikes_binned_ds,ds_rate);
% this will be inputs to the function
n = 1; % neuron #
% not sure where to put this
vxN = [vxN;0];
vyN = [vyN;0];
r = [r;0];
phi = [phi;0];

%% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:max(spikes_binned)
    ind = [ind; find(spikes_binned(:,n) >= numspikes)];
end;

figure('Name','Occupancy Normalized Histogram');
warning('off');

%% x-velocity
% Histogram of spiking to x-velocity
% subplot(2,2,1);
    velocities = min(vxN):range(vxN)/80:max(vxN);
%     velocities = -.04:.001:.04;
    bar(velocities,hist(vxN(ind),velocities)./hist(vxN,velocities));
    xlabel('x velocity');
    ylabel('normalized spike counts');

    % GLM
    cov1 = [vxN];
    b1 = glmfit(cov1,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b1(1)+b1(2)*velocities),'r');
    
    
    % model 2
    cov2 = [vxN.^2];
    b2 = glmfit(cov2,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b2(1)+b2(2)*velocities.^2),'g');
    
    % model 3
    cov3 = [vxN vxN.^2];
    b3 = glmfit(cov3,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b3(1)+b3(2)*velocities+b3(3)*velocities.^2),'c');
    
    % model 4
    cov4 = [vxN.^3];
    b4 = glmfit(cov4,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b4(1)+b4(2)*velocities.^3),'m');
    
    % model 5
    cov5 = [vxN vxN.^3];
    b5 = glmfit(cov5,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b5(1)+b5(2)*velocities+b5(3)*velocities.^3),'b');
    
    % model 6
    cov6 = [vxN vxN.^2 vxN.^3];
    b6 = glmfit(cov6,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b6(1)+b6(2)*velocities+b6(3)*velocities.^2+b6(4)*velocities.^3),'y');
    
    % legend
    legend('velocities','lin','quad','lin+quad','cub','li+cub','lin+quad+cub')
    set(gca, 'fontsize',14)
    
%     % error calculation(?)
%     lambdaEst{1} = gen_lambda(b1,cov1);
%     lambdaEst{2} = gen_lambda(b2,cov2);
%     lambdaEst{3} = gen_lambda(b3,cov3);
%     lambdaEst{4} = gen_lambda(b4,cov4);
%     lambdaEst{5} = gen_lambda(b5,cov5);
%     lambdaEst{6} = gen_lambda(b6,cov6);
%     
%     plot_ks(spikes_binned,lambdaEst);
%     cur_title = get(gca, 'Title');
%     title([cur_title.String ': neuron ' num2str(n)]);
% %     saveas(gcf, ['KS-neuron_' num2str(n) '.png'])

    
%% y-velocity
% Histogram of spiking to y-velocity 
subplot(2,2,2);
    velocities = min(vyN):range(vyN)/80:max(vyN);
    bar(velocities,hist(vyN(ind),velocities)./hist(vyN,velocities));
    xlabel('y velocity');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(vyN,spikes_binned(:,n),'poisson');
    hold on;
    plot(velocities,exp(b(1)+b(2)*velocities),'r');

%% movement speed
% Histogram of spiking to movement speed 
subplot(2,2,3);
    rs = 0:range(r)/40:max(r);
%     rs = 0:.001:.04;
    bar(rs,hist(r(ind),rs)./hist(r,rs));
    xlabel('speed');
    ylabel('normalized spike counts');
    
    % GLM
    b = glmfit(r,spikes_binned(:,n),'poisson');    
    hold on;
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
    plot(phis,exp(b(1)+b(2)*phis),'r');
