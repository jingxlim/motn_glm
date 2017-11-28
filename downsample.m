% downsample.m
% This function downsamples the data to reduce computation time
% rate is desired sampling rate in Hz. It must produce an integer timestep
% in milliseconds

function [xDS,yDS,spikesDS] = downsample(x,y,spikes,rate)

step = 1/rate*1000;

if rem(step,1)
    error('Invalid sampling rate. Please choose a sampling rate that produces an integer timestep in ms.');
end

xDS = x(1:step:end);
yDS = y(1:step:end);

spikesDS = NaN(length(xDS),size(spikes,2));
for n = 1:size(spikes,2)
    for i = 1:ceil(length(x)/step)
        if i < ceil(length(x)/step)
            spikesDS(i,n) = sum(spikes(1+(i-1)*step:i*step,n));
        else
            spikesDS(i,n) = sum(spikes(1+(i-1)*step:end,n));
        end
    end
end

end