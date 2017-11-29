function [ISI,XLocAtSpikes,YLocAtSpikes] = gen_spike_stat(x,y,spikes)

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