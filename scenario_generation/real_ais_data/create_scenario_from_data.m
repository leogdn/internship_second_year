function [Y,Z,Y_cart,Z_cart,xMarkov,trajWithSpoofing,duration] = create_scenario_from_data(rawData,F,probSpoofedTraj,Sigma,durSpoofMin,angle_comp,theta,posAntenna,radius)

K = size(Sigma,1);

[indsChosen,duration] = choose_duration(rawData,"quartile");
dataLong = rawData(indsChosen);
dataCut = cut_series(dataLong,duration);
dT = choose_dT(dataCut);
[tLine,dataInterp] = interpolate_data(dataCut,dT,duration); 
data = preprocess_data(dataInterp,posAntenna,radius,false); %with velocities = [theta,rho,dTheta,dRho] ; without velocities = [theta,rho]

T = length(tLine);
N=length(dataInterp);
Y = cell(1,N);
Z = cell(1,N);
Y_cart = cell(1,N);
Z_cart = cell(1,N);
cellx = cell(1,N);

d = size(data{1},2);

Id = eye(d);

bar = waitbar(0,'Generating trajectories');

xMarkov = cell(1,N);

trajWithSpoofing = cell(1,N);

for n=1:N 
    trajWithSpoofing{n} = rand<probSpoofedTraj;
    xMarkov{n}=ones(1,T);
    xMarkov{n}(1) = 1;
    Y{n} = zeros(T,d);
    Z{n} = zeros(T,d);
    changeXShort = zeros(1,T);
    if trajWithSpoofing{n}
        for t = 2:T
            xMarkov{n}(t) = randsample(1:K, 1, 'true', F(xMarkov{n}(t-1), :));
            changeXShort(t) = (xMarkov{n}(t) ~= xMarkov{n}(t-1));
        end
    end
    changeXShort(T) = 1;
    changeXShort(1) = 1;
    tChangeShort = find(changeXShort);
    nChangeShort = length(tChangeShort);
    changeX = changeXShort;

    for iC=1:nChangeShort-1
        tC = tChangeShort(iC);
        tNextC = tChangeShort(iC+1);
        if (tNextC - tC < durSpoofMin) && (xMarkov{n}(tC) == 2)
            changeX(tC) = 0;
            changeX(tNextC) = 0;
            xMarkov{n}(tC:tNextC) = ones(1,tNextC - tC + 1);
        end
    end
    for t = 1:T
        % Y{n}(t,:) = data{n}(t,:) + sqrt(SigmaK(xMarkov(t), :)).*randn(1,d);
        % Z{n}(t,:) = data{n}(t,:) + sqrt(SigmaR).*randn(1, d) + theta*Id(angle_comp, :);
        Y{n}(t,:) = data{n}(t,:);
        Z{n}(t,:) = data{n}(t,:) + sqrt(Sigma(xMarkov{n}(t), :)).*randn(1,d) - theta*Id(angle_comp, :);
    end
    %disp(n);

    [YxArray,YyArray] = pol2cart(Y{n}(:,2),Y{n}(:,1));
    Y_cart{n}(:,1) = YxArray;
    Y_cart{n}(:,2) = YyArray;
    [ZxArray,ZyArray] = pol2cart(Z{n}(:,2),Z{n}(:,1));
    Z_cart{n}(:,1) = ZxArray;
    Z_cart{n}(:,2) = ZyArray;

    waitbar(n/N, bar);
end

end