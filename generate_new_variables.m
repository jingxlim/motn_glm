% generate_new_variables.m
% [ISI,XLocAtSpikes,YLocAtSpikes] = generate_new_variables(x,y,spikes)
% This function returns cell arrays of Interspike intervals (ISIs), X
% Locations at spike times, and Y Locations at spike times.

function [Vx,Vy,dir,speed] = generate_new_variables(x,y,spikes,pos_sampling_rate)

% n = size(spikes,2);

% for i = 1:n
Vx = diff(x)*pos_sampling_rate;
Vy = diff(y)*pos_sampling_rate;
dir = atan2((y(2:end)-y(1:end-1)),x(2:end)-x(1:end-1));
speed = (Vx.^2+Vy.^2).^.5;
% end

end

