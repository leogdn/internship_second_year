clear; clc; close all; fc = 0;

T = 4000; % number of time steps T=2000
d = 2; % dimension of the observations
angle_comp = 1; % the component where the angle information is stored
theta = 0.1; % angular bias
N = 40; % number of tracks for AIS and Radar
K = 2; % number of modes

%F = 0.9*eye(K) + (0.1/(K-1))*(ones(K) - eye(K)); % Transition matrix
F = [0.999 0.001;
    0.005 0.995];
probSpoofedTraj = 0.3;
%Sigma = (0.005*(1:K)').^2*rand(1, d); % noise variances
sigmaNoise = (.25*rand(1,4));
%sigmaNoise = (0.1.*rand(1,2)).^2;
%sigmaNoise = zeros(1,4);

%Sigma = (0.5*(1:K)').^2*rand(1, d);

Y0 = [1, 1, 1, 1]';
YT = [5,1,5,1]';
transMat = F;
durSpoofMin = 0.1*T;

posAntenna = [];
radius = Inf;
source_dir = "C:\Users\lgaudin\Desktop\ushant_ais\data";
%rawData = get_raw_data(source_dir,100);

%[Y,Z] = create_scenario_random(T,N,d,angle_comp,theta,Sigma,F,K);
[Y,Z,Y_cart,Z_cart,sigmaSpoofing,transMatEst,xMarkov,trajWithSpoofing] = create_scenario_sampling(N,T,Y0,YT,"normal",0,transMat,probSpoofedTraj,durSpoofMin,theta,angle_comp,sigmaNoise);
%[Y,Z,xMarkov,trajWithSpoofing,duration] = create_scenario_from_data(rawData,F,probSpoofedTraj,Sigma,durSpoofMin,angle_comp,theta,posAntenna,radius);

N = size(Y,2);

realAssoc = (1:N)';

T = size(Y{1},1);
  
% initial estimates
Sigma0 = sort(rand(K, 1))*ones(1, d);
theta0 = 0.2;
F = F./sum(F, 2);

F0 = generate_intial_F(K);

diagSigmaSpoofing = diag(sigmaSpoofing);

Sigma = [sigmaNoise([1 3]);diagSigmaSpoofing([1 3])'];

% Sigma = zeros(2,2);
% Sigma(1,:) = sigmaNoise([1 3]);
% Sigma(2,1) = sigmaSpoofing(1,1);
% Sigma(2,2) = sigmaSpoofing(3,3);

B = 5; % the number of best assignments to be considered
a_EM = 0.65; % power for the coefficient of stochastic approximation
%a_EM = .95;
burnin_EM = 0.01*T;
tic;
[Thetas, Fs, Sigmas, meanThetas, meanFs, meanSigmas, errorMats, scoreTables, scoreMeans, FiltX] = T2TA(Z, Y, K, F0, Sigma0, theta0, angle_comp, B, a_EM, burnin_EM,realAssoc); %Estimation of the parameters
toc;

estX = cell(1,N);
filtXMean = cell(1,N);

for n=1:N
    estX{n} = 1*(FiltX{n} > 0.75) + 2*(FiltX{n} <= 0.75);
    %[estX{n},filtXMean{n}] = detect_spoofing(FiltX{n}(:,2),0.75,20);
end

%% Compute the confusion matrix and the F1-score

% Compute the confusion matrix
confusionMatrix = zeros(K,K,T);

for t = 1:T
    for n = 1:N
        trueLabel = xMarkov{n}(t);
        estLabel = estX{n}(t);
        confusionMatrix(trueLabel, estLabel, t) = confusionMatrix(trueLabel, estLabel, t) + 1;
    end
end

%Compute the F1-score
precision = confusionMatrix(1,1,:)./(confusionMatrix(1,1,:)+confusionMatrix(2,1,:));
recall = confusionMatrix(1,1,:)./(confusionMatrix(1,1,:)+confusionMatrix(1,2,:));
f1Score = 2 * (precision .* recall)./(precision + recall);

%% Recompute the means

thresholdMean = 0.25*T;

for t=1:T
    if t==1
        meanThetas(t) = 0;
        meanSigmas(:,:,t) = 0;
        meanFs(:,:,t) = 0;
    else
        meanThetas(t) = update_mean(meanThetas(t-1),Thetas(t),thresholdMean,t);
        meanSigmas(:,:,t) = update_mean(meanSigmas(:,:,t-1),Sigmas(t),thresholdMean,t);
        meanFs(:,:,t) = update_mean(meanFs(:,:,t-1),Fs(t),thresholdMean,t);
    end
end

%% plot the results

refValueLine.color = 'r';
refValueLine.size = 3;

estValueLine.color = 'b';
estValueLine.size = 3;

% figure;
% plot(Thetas,'LineWidth',4);
% hold on;
% plot(theta*ones(1, T), 'r','LineWidth',3);
% hold off;
% 
% fc = fc + 1; figure(fc);
% plot(Thetas);
% hold on;
% plot(theta*ones(1, T), 'r');
% hold off;
% title(sprintf('$\\theta$'), 'interpreter', 'Latex');
% 
% fc = fc + 1; figure(fc);
% for i = 1:K
%     for j = 1:d
%         subplot(K, d, (i-1)*d + j)
%         plot(squeeze(Sigmas(i, j, :)));
%         hold on;
%         plot(Sigma(i, j)*ones(1, T), 'r');
%         hold off;
%         title(sprintf('$\\Sigma(%d, %d)$', i, j), 'interpreter', 'Latex');
%         xlabel('$t$', 'Interpreter', 'Latex');
%     end
% end
% 
% fc = fc + 1; figure(fc);
% for i = 1:K
%     for j = 1:K
%         subplot(K, K, (i-1)*K + j);
%         plot(squeeze(Fs(i, j, :)));
%         hold on;
%         plot(F(i, j)*ones(1, T), 'r');
%         hold off;
%         title(sprintf('$F(%d, %d)$', i, j), 'interpreter', 'Latex')
%     end
% end

numFigs = 7;
fc=0;
figTab = gobjects(1,numFigs);
fc = fc + 1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.2, 11]);
plot(Thetas,'LineWidth',estValueLine.size,'DisplayName','$\theta$');
hold on;
plot(meanThetas,'g','LineWidth',estValueLine.size,'DisplayName','$\theta_M$');
plot(theta*ones(1, T), 'r','LineWidth',refValueLine.size,'DisplayName','Real value');
legend('show','Interpreter','latex');
hold off;
%title(sprintf('$\\theta$'), 'interpreter', 'Latex');
xlabel('$t$', 'Interpreter', 'Latex');
ylabel(sprintf('$\\theta$'), 'interpreter', 'Latex');
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

fc = fc + 1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 34.4, 22]);
for i = 1:K
    for j = 1:d
        subplot(K, d, (i-1)*d + j)
        plot(squeeze(Sigmas(i, j, :)),'LineWidth',estValueLine.size,'DisplayName','$\Sigma$');
        hold on;
        plot(squeeze(meanSigmas(i,j,:)),'g','LineWidth',estValueLine.size,'DisplayName','$\Sigma_M$');
        plot(Sigma(i, j)*ones(1, T), 'r','LineWidth',refValueLine.size,'DisplayName','Real value');
        legend('Interpreter','latex');
        hold off;
        %title(sprintf('$\\Sigma(%d, %d)$', i, j), 'interpreter', 'Latex');
        xlabel('$t$', 'Interpreter', 'Latex');
        ylabel(sprintf('$\\Sigma(%d, %d)$', i, j), 'interpreter', 'Latex');
    end
end
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

fc = fc + 1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 34.4, 22]);
for i = 1:K
    for j = 1:K
        subplot(K, K, (i-1)*K + j);
        plot(squeeze(Fs(i, j, :)),'LineWidth',estValueLine.size,'DisplayName','$F$');
        hold on;
        plot(squeeze(meanFs(i,j,:)),'g','LineWidth',estValueLine.size,'DisplayName','$F_M$');
        plot(F(i, j)*ones(1, T), 'r','LineWidth',refValueLine.size,'DisplayName','Real value');
        legend('Interpreter','latex');
        hold off;
        %title(sprintf('$F(%d, %d)$', i, j), 'interpreter', 'Latex');
        xlabel('$t$', 'Interpreter', 'Latex');
        ylabel(sprintf('$F(%d, %d)$', i, j), 'interpreter', 'Latex');
    end
end
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

% fc = fc + 1; figure(fc);
% plot(meanThetas,'LineWidth',estValueLine.size);
% hold on;
% plot(theta*ones(1, T), 'r','LineWidth',refValueLine.size);
% hold off;
% %title(sprintf('$\widehat{\theta}_M$'), 'interpreter', 'Latex');
% xlabel('$t$', 'Interpreter', 'Latex');
% ylabel(sprintf('$\\widehat{\\theta}_M$'), 'interpreter', 'Latex');
% allAxes = findall(figure(fc), 'Type', 'axes');
% for k = 1:length(allAxes)
%     allAxes(k).FontSize = 14;
% end
% figTab(fc) = gcf;

% fc = fc + 1; figure(fc);
% for i = 1:K
%     for j = 1:d
%         subplot(K, d, (i-1)*d + j);
%         plot(squeeze(meanSigmas(i, j, :)),'LineWidth',estValueLine.size);
%         hold on;
%         plot(Sigma(i, j)*ones(1, T), 'r','LineWidth',refValueLine.size);
%         hold off;
%         %title(sprintf('$\\Sigma_M(%d, %d)$', i, j), 'interpreter', 'Latex');
%         xlabel('$t$', 'Interpreter', 'Latex');
%         ylabel(sprintf('$\\Sigma_M(%d, %d)$', i, j), 'interpreter', 'Latex');
%     end
% end
% allAxes = findall(figure(fc), 'Type', 'axes');
% for k = 1:length(allAxes)
%     allAxes(k).FontSize = 15;
% end
% figTab(fc) = gcf;

% fc = fc + 1; figure(fc);
% for i = 1:K
%     for j = 1:K
%         subplot(K, K, (i-1)*K + j);
%         plot(squeeze(meanFs(i, j, :)),'LineWidth',estValueLine.size);
%         hold on;
%         plot(F(i, j)*ones(1, T), 'r','LineWidth',refValueLine.size);
%         hold off;
%         %title(sprintf('$F_M(%d, %d)$', i, j), 'interpreter', 'Latex');
%         xlabel('$t$', 'Interpreter', 'Latex');
%         ylabel(sprintf('$F_M(%d, %d)$', i, j), 'interpreter', 'Latex');
%     end
% end
% allAxes = findall(figure(fc), 'Type', 'axes');
% for k = 1:length(allAxes)
%     allAxes(k).FontSize = 15;
% end
% figTab(fc) = gcf;

fc = fc + 1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.2, 11]);
plot(squeeze(scoreTables(1,:)),'LineWidth',estValueLine.size);
xlabel('$t$', 'Interpreter', 'Latex');
ylabel("Error");
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

fc = fc + 1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.2, 11]);
plot(squeeze(scoreMeans(:,1)),'LineWidth',estValueLine.size);
xlabel('$t$', 'Interpreter', 'Latex');
ylabel("Error");
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

% fc = fc + 1; figure(fc);
% for i = 1:B
%     subplot(1,B,i);
%     plot(squeeze(scoreTables(i,:)),'LineWidth',estValueLine.size);
%     xlabel('$t$', 'Interpreter', 'Latex');
%     ylabel("Error");
% end
% allAxes = findall(figure(fc), 'Type', 'axes');
% for k = 1:length(allAxes)
%     allAxes(k).FontSize = 15;
% end
% figTab(fc) = gcf;

% fc = fc + 1; figure(fc);
% for i = 1:B
%     plot(scoreTables(i,:),'LineWidth',estValueLine.size);
%     xlabel('$t$', 'Interpreter', 'Latex');
%     ylabel("Averaged Error");
% end
% allAxes = findall(figure(fc), 'Type', 'axes');
% for k = 1:length(allAxes)
%     allAxes(k).FontSize = 15;
% end
% figTab(fc) = gcf;

% spoofedTrajTab = find(cell2mat(trajWithSpoofing));
% 
% tracksTab = spoofedTrajTab;
% for nTrack = 1:length(tracksTab)
%     fc = fc + 1; figure(fc);
%     set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.2, 11]);
%     %subplot(3,2,1);
%     track = tracksTab(nTrack);
%     %plot(estX{track}(:),'LineWidth',estValueLine.size);
%     plot(FiltX{track}(:,2),'LineWidth',estValueLine.size,'DisplayName','$\alpha_t(2)$');
%     hold on;
%     %plot(filtXMean{track},'LineWidth',estValueLine.size,'DisplayName','Averaged $\alpha_t(2)$')
%     plot(xMarkov{track}(:),'LineWidth',refValueLine.size,'Color',refValueLine.color,'DisplayName','$x_t$');
%     plot(0.75.*ones(1,T),'m','LineWidth',refValueLine.size,'DisplayName','Threshold');
%     xlabel('$t$', 'Interpreter', 'Latex');
%     legend('Interpreter','latex');
%     hold off;
%     % subplot(3,2,2);
%     % track = tracksTab(nTrack);
%     % %plot(estX{track}(:),'LineWidth',estValueLine.size);
%     % plot(FiltX{track}(:,2),'LineWidth',estValueLine.size,'DisplayName','$\alpha_t(2)$');
%     % hold on;
%     % plot(xMarkov{track}(:),'LineWidth',refValueLine.size,'Color',refValueLine.color,'DisplayName','$x_t$');
%     % xlabel('$t$', 'Interpreter', 'Latex');
%     % legend('Interpreter','latex');
%     % hold off;
%     % subplot(3,2,3);
%     % plot(Y_cart{track}(1,:),'LineWidth',estValueLine.size,'DisplayName',"AIS");
%     % hold on;
%     % plot(Z_cart{track}(1,:),'LineWidth',estValueLine.size,'DisplayName',"Radar");
%     % ylabel('$x$','Interpreter', 'Latex');
%     % xlabel('$t$', 'Interpreter', 'Latex');
%     % legend;
%     % hold off;
%     % subplot(3,2,4);
%     % plot(Y{track}(:,angle_comp),'LineWidth',estValueLine.size,'DisplayName',"AIS");
%     % hold on;
%     % plot(Z{track}(:,angle_comp),'LineWidth',estValueLine.size,'DisplayName',"Radar");
%     % ylabel('$\theta$','Interpreter', 'Latex');
%     % xlabel('$t$', 'Interpreter', 'Latex');
%     % legend;
%     % hold off;
%     % subplot(3,2,5);
%     % plot(Y_cart{track}(3,:),'LineWidth',estValueLine.size,'DisplayName',"AIS");
%     % hold on;
%     % plot(Z_cart{track}(3,:),'LineWidth',estValueLine.size,'DisplayName',"Radar");
%     % ylabel('$y$','Interpreter', 'Latex');
%     % xlabel('$t$', 'Interpreter', 'Latex');
%     % legend;
%     % hold off;
%     % subplot(3,2,6);
%     % plot(Y{track}(:,2),'LineWidth',estValueLine.size,'DisplayName',"AIS");
%     % hold on;
%     % plot(Z{track}(:,2),'LineWidth',estValueLine.size,'DisplayName',"Radar");
%     % ylabel('$r$','Interpreter', 'Latex');
%     % xlabel('$t$', 'Interpreter', 'Latex');
%     % legend;
%     % hold off;
% 
%     allAxes = findall(figure(fc), 'Type', 'axes');
%     for k = 1:length(allAxes)
%         allAxes(k).FontSize = 15;
%     end
%     figTab(fc) = gcf;
% end

fc = fc+1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 34.4, 22]);
for i=1:size(confusionMatrix,1)
    for j=1:size(confusionMatrix,2)
        subplot(2,2,(i-1)*2+j);
        plot(squeeze(confusionMatrix(i,j,:)),'LineWidth',estValueLine.size);
        title(sprintf('(%d,%d)',i,j));
        disp(i);
        disp(j);
    end
end
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

fc = fc+1; figure(fc);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 17.2, 11]);
plot(squeeze(f1Score(:)),'LineWidth',estValueLine.size);
title(sprintf('F1-score'));
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end
figTab(fc) = gcf;

export_figures(figTab,'figures',T);