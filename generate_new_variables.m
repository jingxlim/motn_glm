% generate_new_variables.m
% [ISI,XLocAtSpikes,YLocAtSpikes] = generate_new_variables(x,y,spikes)
% This function returns cell arrays of Interspike intervals (ISIs), X
% Locations at spike times, Y Locations at spike times, movement direction,
% and movement speed.

function [Vx,Vy,dir,speed] = generate_new_variables(x,y,pos_sampling_rate)

% add zero for length
Vx    = [0; diff(x)*pos_sampling_rate];
Vy    = [0; diff(y)*pos_sampling_rate];
speed = (Vx.^2+Vy.^2).^.5;

% repeat first index instead of zero
dir   = atan2((y(2:end)-y(1:end-1)),x(2:end)-x(1:end-1));
dir   = [dir(1); dir];

end

