function dT = choose_dT(data)
    minArray = cellfun(@minPeriod,data);
    dTMin = min(minArray);
    dT = (1/2)*dTMin; % Assign the minimum time difference to dT
end

function minT = minPeriod(mat)
    tLine = mat(:,5);
    dtLine = diff(tLine);
    minT = min(dtLine);
end