function [Y,Z] = create_scenario_sampling2(N,T,Y0,YT,transMat,theta,compAngle)
    
    x_vec = linspace(-3, 3, 10);
    
    t_given = linspace(1, T, 10);
    
    X_given = num2cell(10*normpdf(x_vec, 0, 1));
    % X_given = {0, 10, 20, 30, 30, 20, 0};
    dim_given = repmat({1}, 1, 10);
    
    %mu_0 = [0, 1, 0, 1]';
    %Cov_0 = 0.1*eye(4);
    mu_0 = Y0;
    Cov_0 = 0.1*eye(4);

    F = [[1, 0.1, 0, 0]; [0, 0.99, 0, 0]; [0, 0, 1, 0.1]; [0, 0, 0, 0.99]];
    U = diag([0.001, 0.01, 0.001, 0.01]);
    
    Y_cart = cell(1,N);
    Z_cart = cell(1,N);
    Y = cell(1,N);
    Z = cell(1,N);

    K=2;

    for n = 1:N
        %%Modeling of spoofing

        giveTSpoofed(1) = 1;
        dim_given_spoofed{1} = [1, 2, 3, 4];
        X_given_spoofed{1} = X_genuine(:,1)';
        giveTSpoofed(T) = 1;
        dim_given_spoofed{T} = [1, 2, 3, 4];
        X_given_spoofed{T} = X_genuine(:,T)';
        t_given_spoofed = find(giveTSpoofed);
        dim_given_spoofed = dim_given_spoofed(~cellfun(@isempty,dim_given_spoofed));
        X_given_spoofed = X_given_spoofed(~cellfun(@isempty,X_given_spoofed));

        X_spoofed = conditional_sampling(mu_0, Cov_0, F, U, T, t_given_spoofed, X_given_spoofed, dim_given_spoofed);

        Y_cart{n} = X_genuine;
        Z_cart{n} = X_spoofed;

        Y{n} = Y_cart{n};
        Z{n} = Z_cart{n};

        I_d = eye(4);

        % Y{n} = convert_cartesian_to_polar_sampling(Y_cart{n});
        % Z{n} = convert_cartesian_to_polar_sampling(Z_cart{n});
        % Z{n} = Z{n} + theta*I_d(compAngle,:); %Check if we have to add bias to angular velocity
        
        %+ theta*I_d(:,comp_theta)

        %Uncomment to get only the position 
        % Y{n}(:,2) = [];
        % Y{n}(:,3) = [];
        % Z{n}(:,2) = [];
        % Z{n}(:,3) = [];

    end
    
    %X0 = conditional_sampling(mu_0, Cov_0, F, U, T, [], {}, {});
    
    
    % t_given = [1 T];
    % X_given = {[X0(1, 1) -X0(2, 1) X0(3, 1) -X0(4, 1)], [X0(1, T) -X0(2, T) X0(3, T) -X0(4, T)]'};
    % dim_given = {[1, 2, 3, 4], [1, 2, 3, 4]};
    % 
    % X = conditional_sampling(mu_0, Cov_0, F, U, T, t_given, X_given, dim_given);
    % Y = conditional_sampling(mu_0, Cov_0, F, U, T, [T], {0}, {1});
    % scatter(X0(1, :), X0(3, :), 'b');
    % hold on;
    % scatter(X(1, :), X(3, :), 'r');
    % xlabel('X');
    % ylabel('Y');
    
    % T_medium = floor(T/2);
    % 
    % t_given_genuine = [1 T_medium T];
    % dim_given_genuine = {[1, 2, 3, 4], [1, 3], [1, 2, 3, 4]};
    % X_middle_genuine = 0.5.*([X0(1, 1) X0(3, 1)]+[X0(1, T) X0(3, T)]);
    % X_given_genuine = {[X0(1, 1) -X0(2, 1) X0(3, 1) -X0(4, 1)], X_middle_genuine, [X0(1, T) -X0(2, T) X0(3, T) -X0(4, T)]'};
    % 
    % t_given_spoofed = [1 T_medium T];
    % dim_given_spoofed = {[1, 2, 3, 4], [1, 3], [1, 2, 3, 4]};
    % [X_middle(1),X_middle(2)] = point_mediator([X0(1,1);X0(3,1)],[X0(1,T);X0(3,T)],0.5);
    % X_given_spoofed = {[X0(1, 1) -X0(2, 1) X0(3, 1) -X0(4, 1)], [X_middle(1) X_middle(2)], [X0(1, T) -X0(2, T) X0(3, T) -X0(4, T)]'};
    % 
    % X_genuine = conditional_sampling(mu_0, Cov_0, F, U, T, t_given_genuine, X_given_genuine, dim_given_genuine);
    % X_spoofed = conditional_sampling(mu_0, Cov_0, F, U, T, t_given_spoofed, X_given_spoofed, dim_given_spoofed);
    
    %scatter(X0(1, :), X0(3, :), 'b');
    % hold on;
    % scatter(X_middle(1),X_middle(2),'k');
    % scatter(X_genuine(1, :), X_genuine(3, :), 'g');
    % scatter(X_spoofed(1, :), X_spoofed(3, :), 'r');
    % plot(X_middle_genuine(1),X_middle_genuine(2),'bo', 'MarkerSize', 10);
    % plot(X_middle(1),X_middle(2),'b*', 'MarkerSize', 10)
    % xlabel('X');
    % ylabel('Y');
    % hold off;
    
    %%
    
    % delta_t_vec = 0.1*ones(1, T);
    % v_vec = X0([2, 4], :);
    % gamma_vec = [0.9, 0.9];
    % Sigma = diag([0.1, 0.1]);
    % Theta = diag(gamma_vec);
    % X_OU = OU_process(delta_t_vec, v_vec, gamma_vec, Theta, Sigma);
    % 
    % scatter(X_OU(1, :), X_OU(2, :), 'k');
    % hold off;

end