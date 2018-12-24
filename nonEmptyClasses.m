function n = nonEmptyClasses(data, cutoff)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    classes = unique(data);
    dataMatrix = repmat(data,[1 1 numel(classes)]);
    classMatrix = permute(repmat(classes,[1,size(data,1),size(data,2)]),[2 3 1]);
    counts = squeeze((sum(dataMatrix == classMatrix,2)./size(data,2))>cutoff);
    n=mode(sum(counts,2));
end

