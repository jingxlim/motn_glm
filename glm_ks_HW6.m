%% Models of the Neuron: Evaluating a statistical model (part 2)
% Homework 6
% Script glm_part1_ks.m
%
% MATLAB code to fit a GLM model of the relation between
% spiking and the rat's position, and draw K-S plot for the
% second Neuroinformatics 2005 GLM problem set.

% load the rat trajectory and spiking data;
clr;
load('train.mat');
[ISI,x_at_spiketimes,y_at_spiketimes,Vx,Vy,dir] = generate_new_variables(xN,yN,spikes_binned);

% pick a neuron
spikes_binned = spikes_binned(:,8);

% visualize the raw data
figure(7); clf; hold on;
title('Visualize the raw data');
plot(xN,yN,x_at_spiketimes{8},y_at_spiketimes{8},'r.');

%% 4. Characterize the goodness-of-fit between data and model
% fit a GLM model to the x and y positions. (ADD ALL YOUR MODEL CANDIDATES
% HERE!!! Two possible models are listed below.

[b_lin,dev_lin,stats_lin] = glmfit([xN yN xN.^2 yN.^2],spikes_binned,'poisson');
[b_quad,dev_quad,stats_quad] = glmfit([xN.^2 yN.^2],spikes_binned,'poisson');

%*******  K-S Plot  *******************
% Note that the ks plot for one model is given below. You should overlay ks
% plots for each of your models.

% graph the K-S plot and confidence intervals for the K-S statistic

% first generate the conditional intensity at each timestep
% ** Adjust the below line according to your choice of model.
% remember to include a column of ones to multiply the default constant GLM
% parameter beta_0**

 % based on our GLM model with the log "link function"
                       % Use your parameter estimates (b) from glmfit along
                       % with the covariates you used (xN, yN, ...)

lambdaEst_lin = exp(b_lin(1)+b_lin(2)*xN+ b_lin(3)*yN+b_lin(4)*xN.^2+ b_lin(5)*yN.^2);
lambdaEst_quad = exp(b_quad(1)+b_quad(2)*xN.^2+ b_quad(3)*yN.^2);

% linear-quadratic model
timestep = 1;
lambdaInt = 0;
j=0;

for t=1:length(spikes_binned)
    lambdaInt = lambdaInt + lambdaEst_lin(t)*timestep;
    if (spikes_binned(t))
        j = j + 1;
        KS_lin(j) = 1-exp(-lambdaInt);
        lambdaInt = 0;
    end
end

KSSorted_lin = sort(KS_lin);
N_lin = length(KSSorted_lin);

% quadratic model
timestep = 1;
lambdaInt = 0;
j=0;

for t=1:length(spikes_binned)
    lambdaInt = lambdaInt + lambdaEst_quad(t)*timestep;
    if (spikes_binned(t))
        j = j + 1;
        KS_quad(j) = 1-exp(-lambdaInt);
        lambdaInt = 0;
    end
end

KSSorted_quad = sort(KS_quad);
N_quad = length(KSSorted_quad);

% plot "perfect" models
figure(8); clf; hold on;
plot(0:.01:1,0:.01:1, 'g',...
    0:.01:1, [0:.01:1]+1.36/sqrt(N_lin),'r',...
    0:.01:1,[0:.01:1]-1.36/sqrt(N_lin), 'r')

% plot KS plots
h_lin = plot(([1:N_lin]-.5)/N_lin, KSSorted_lin, 'b',...
    'DisplayName', 'Linear-quadratic model');
h_quad = plot(([1:N_quad]-.5)/N_quad, KSSorted_quad, 'k',...
    'DisplayName', 'Quadratic model');

axis( [0 1 0 1] );
xlabel('Uniform CDF');
ylabel('Empirical CDF of Rescaled ISIs');
title('KS Plot with 95% Confidence Intervals');
legend([h_lin h_quad], 'Location', 'northwest');

ks_stat_lin = max(abs(KSSorted_lin - ([1:N_lin]-.5)/N_lin))
ks_stat_quad = max(abs(KSSorted_quad - ([1:N_quad]-.5)/N_quad))

%%
% The linear-quadratic model does a better job than the purely quadratic
% model in capturing the statistical structure of the data, since it has a 
% lower K-S statistic (0.1132 vs 0.2357), which represents the maximum
% difference between the cumulative distribution of the model to an
% independent identically distributed exponential set of random variables.
%
% That said, the linear-quadratic model still does not capture all of the
% statistical structure of the data, though it is very close. This can be
% seen from the K-S plot, where the distribution of the linear-quadratic
% model deviates significantly from the unity (45 degree) line (which
% represents perfect fit of the model to the data), at multiple locations.
% Only the regions that sits within the 95 confidence bounds can be
% considered good fits.