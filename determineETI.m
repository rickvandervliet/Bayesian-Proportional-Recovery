function output = determineETI(samples, alpha)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

samples = samples(~isnan(samples));
output = NaN(3,1);

if ~isempty(samples)
    sortedSamples = sort(samples);
    interval(1) = sortedSamples(max([1 round(numel(samples)*alpha/2)]));
    interval(2) = sortedSamples(round(numel(samples)*(1-alpha/2)));
    output = [mean(samples) interval];
end;    