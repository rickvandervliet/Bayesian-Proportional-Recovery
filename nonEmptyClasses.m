function [n,p] = nonEmptyClasses(data, cutoff)
    n = median(sum(data>cutoff,2));

    p = NaN(size(data,1),1);
    for ii=1:size(data,1)
        p(ii) = sum(data(ii,data(ii,:)>cutoff));
    end;
    p = mean(p);
end
