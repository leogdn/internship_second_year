function [path_genuine,path_spoofed] = create_path_sampling(T,Y0,YT)
    
    x_vec = linspace(-3, 3, 10);
    
    t_given = linspace(1, T, 10);
    
    %X_given = num2cell(10*normpdf(x_vec, 0, 1));
    % X_given = {0, 10, 20, 30, 30, 20, 0};
    dim_given = repmat({1}, 1, 10);
    
    %mu_0 = [0, 1, 0, 1]';
    %Cov_0 = 0.1*eye(4);
    mu_0 = Y0;
    Cov_0 = 0.1*eye(4);

    F = [[1, 0.1, 0, 0]; [0, 0.99, 0, 0]; [0, 0, 1, 0.1]; [0, 0, 0, 0.99]];
    U = diag([0.001, 0.01, 0.001, 0.01]);
    

    T_medium = floor(T/2);
    t_given_genuine = [1 T_medium T];
    dim_given_genuine = {[1, 2, 3, 4], [1, 3], [1, 2, 3, 4]};
    X_middle_genuine = 0.5.*([Y0(1) Y0(3)]+[YT(1) YT(3)]);
    X_given_genuine = {[Y0(1) Y0(2) Y0(3) Y0(4)], X_middle_genuine, [YT(1) YT(2) YT(3) YT(4)]};
    
    t_given_spoofed = [1 T_medium T];
    dim_given_spoofed = {[1, 2, 3, 4], [1, 3], [1, 2, 3, 4]};
    X_middle_spoofed = zeros(1,2);
    StartPoint = [Y0(1) Y0(3)];
    EndPoint = [YT(1) YT(3)];
    X_middle_spoofed = point_mediator(StartPoint,EndPoint,0.7);
    X_given_spoofed = {[Y0(1) Y0(2) Y0(3) Y0(4)], X_middle_spoofed, [YT(1) YT(2) YT(3) YT(4)]};
    
    path_genuine = conditional_sampling(mu_0, Cov_0, F, U, T, t_given_genuine, X_given_genuine, dim_given_genuine);
    path_spoofed = conditional_sampling(mu_0, Cov_0, F, U, T, t_given_spoofed, X_given_spoofed, dim_given_spoofed);

    %X0 = conditional_sampling(mu_0, Cov_0, F, U, T, [], {}, {});

end