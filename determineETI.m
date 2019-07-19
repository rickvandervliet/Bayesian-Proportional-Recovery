function output = determineETI(samples, alpha)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if isrow(samples)
    samples = samples';
end;    
output = NaN(size(samples,2),3);

for ii=1:size(samples,2)
    sampleset = samples(~isnan(samples(:,ii)),ii);

    if ~isempty(sampleset)
        interval(1) = prctile(sampleset,100*(alpha/2));
        interval(2) = prctile(sampleset,100*(1-alpha/2));
        output(ii,:) = [mean(sampleset) interval];
    end;
end;    