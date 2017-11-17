% generate_new_variables.m
% [ISI,XLocAtSpikes,YLocAtSpikes] = generate_new_variables(x,y,spikes)
% This function returns cell arrays of Interspike intervals (ISIs), X
% Locations at spike times, and Y Locations at spike times.

function [ISI,XLocAtSpikes,YLocAtSpikes] = generate_new_variables(x,y,spikes)

n = size(spikes,2);

ISI = cell(n,1);
XLocAtSpikes = cell(n,1);
YLocAtSpikes = cell(n,1);

for i = 1:n
    ISI{i} = diff(find(spikes(:,i)));
    XLocAtSpikes{i} = x(spikes(:,i));
    YLocAtSpikes{i} = y(spikes(:,i));
end

end

