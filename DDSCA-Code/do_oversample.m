function [X, Y, Z] = do_oversample(X, Y, label)
% Oversampling for class balancing
%
% Date: 2020-01-17
% Author: Lei Du

[n, n_class] = size(label);

for i = 1:n_class
    idxt = find(label(:,i)==1);
    idxf = label(:,i)~=1;
    idx_pool = repmat(idxt,15,1); % repeat the small class for 5 times
    n_diff = n-length(idxt)*1;
    idx_oversample = idx_pool(1:n_diff); % select n_diff subjects for oversampling
    tempX = [X(idxt,:); X(idx_oversample,:); X(idxf,:)];
    tempY = [Y(idxt,:); Y(idx_oversample,:); Y(idxf,:)];
    tempZ = [label(idxt,i); label(idx_oversample,i); label(idxf,i)];
end

X = tempX; Y = tempY; Z = tempZ;