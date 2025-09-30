function [tLine,dataInterp] = interpolate_data(rawData,dT,duration)
    tLine = 0:dT:duration;
    N = length(rawData);
    dataInterp = cellfun(@(mat) interpolate_matrix(mat,tLine),rawData, 'UniformOutput',false);
end

function matInterp = interpolate_matrix(mat,tLine)
    T = size(tLine,2);
    M = size(mat,2);
    matInterp = zeros(T, M); % Initialize the interpolated matrix
    for i = 1:4
        matInterp(:, i+1) = interp1(mat(:, 5), mat(:, i), tLine, 'linear'); % Interpolate each column %%ERROR to check
    end
    matInterp(:,1) = tLine';
end