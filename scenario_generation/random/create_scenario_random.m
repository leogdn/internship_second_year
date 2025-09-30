function [Y,Z] = create_scenario_random(T,N,d,angle_comp,theta,Sigma,F,K)
    Id = eye(d); 
    % simulate the scenarios
    Z = cell(1, N);
    Y = cell(1, N);
    %Each element of these cell is a time serie, the sizes of the cells are (1,N)
    for n = 1:N
        x = zeros(1, T);
        x(1) = 1;    
        for t = 2:T
            x(t) = randsample(1:K, 1, 'true', F(x(t-1), :)); %Selection of the mode of the system, knowing the mode at time t and the transition matrix
        end
        Z{n} = randn*10 + randn(T, d);
        Y{n} = Z{n} + theta*Id(angle_comp, :) + sqrt(Sigma(x, :)).*randn(T, d);
    end
end