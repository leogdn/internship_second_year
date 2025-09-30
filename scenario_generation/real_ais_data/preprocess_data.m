function [data] = preprocess_data(rawData,posAntenna,radius,withVelocities)

% Check if we need posAntenna

n=length(rawData);
data=cell(n,1);

for i=1:n
    %data{i} = rawData{i}(:,components);
    %data{i} = convert_lat_long_to_cartesian(data{i},posAntenna);
    %data{i} = convert_cartesian_to_polar(data{i});
    data{i} = convert_cartesian_to_polar_data(rawData{i});

    indInRange = data{i}(:,2) <= radius;
    %disp(length(indInRange));
    data{i} = data{i}(indInRange,:); %We only keep data that is in range
    if ~withVelocities
        data{i}(:,4) = [];
        data{i}(:,4) = [];
    end
    data{i}(:,1) = [];
end

data = data(cellfun(@isempty, data) == 0);
end

