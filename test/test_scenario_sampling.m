clear all;close all;

Y0 = [1, 1, 1, 1]';
YT = [100,1,100,1]';
N = 60;
K=2;
% transMat = [0.98 0.02;
%             0.01 0.99];
%transMat = 0.9*eye(K) + (0.1/(K-1))*(ones(K) - eye(K));
transMat = [0.999 0.001;
    0.005 0.995];
T = 1000;
%theta = 0.1;
theta = 0.1;
angle_comp = 3;
%sigmaNoise = 0.25.*[rand 0 rand 0];
sigmaNoise = ([ 1 1 .002 .002].*rand(1,4)); %Usual std for SAR radar

probSpoofedTraj = 1.;

durSpoofMin = 0.2*T;

%[Y,Z,Y_cart,Z_cart,sigmaEst,transMatEst,xMarkov] = create_scenario_sampling(N,T,Y0,YT,transMat,probSpoofedTraj,durSpoofMin,theta,comp_theta,sigmaNoise);
[Y,Z,Y_cart,Z_cart,sigmaEst,transMatEst,xMarkov,trajWithSpoofing] = create_scenario_sampling(N,T,Y0,YT,"normal",0,transMat,probSpoofedTraj,durSpoofMin,theta,angle_comp,sigmaNoise);

% for n = 1:N
%     figure;
%     hold on;
%     scatter(Y{n}(1,:),Y{n}(3,:),'r');
%     scatter(Z{n}(1,:),Z{n}(3,:),'g');
%     hold off;
% end
% for n = 1:N
%     figure;
%     hold on;
%     scatter(Y{n}(:,1),Y{n}(:,2),'r');
%     scatter(Z{n}(:,1),Z{n}(:,2),'g');
%     hold off;
% end

refValueLine.color = 'r';
refValueLine.size = 3;

estValueLine.color = 'b';
estValueLine.size = 3;

spoofedTrajTab = find(cell2mat(trajWithSpoofing));
tracksTab = spoofedTrajTab;

fc=0;
fc = fc + 1; figure(fc);

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 34, 22]);

nSamples = 5;
samplesTab = randsample(1:N,nSamples);

for nTrackSamples = 1:nSamples
    hold on;
    nTrack = samplesTab(nTrackSamples);
    if nTrackSamples == 1
        % Afficher les noms pour la légende uniquement au premier passage
        % Stocker les objets graphiques à afficher dans la légende
        hAIS   = scatter(Y_cart{nTrack}(:,1), Y_cart{nTrack}(:,2), 75, 'r', "+");
        hRadar = scatter(Z_cart{nTrack}(:,1), Z_cart{nTrack}(:,2), 75, 'g', "o");
    else
        % Sinon, ne pas définir 'DisplayName'
        scatter(Y_cart{nTrack}(:,1), Y_cart{nTrack}(:,2), 75, 'r', "+");
        scatter(Z_cart{nTrack}(:,1), Z_cart{nTrack}(:,2), 75, 'g', "o");
    end
end

ylabel('$y$', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
lgd = legend([hAIS, hRadar], {'AIS', 'Radar'});
lgd.FontSize = 15;
hold off;

% Met à jour la taille de police pour tous les axes de la figure fc
allAxes = findall(figure(fc), 'Type', 'axes');
for k = 1:length(allAxes)
    allAxes(k).FontSize = 15;
end

% 
% for nTrack = 1:N
%     hold on;
%     scatter(Y_cart{nTrack}(1,:),Y_cart{nTrack}(3,:),75,'r',"+",'DisplayName',"AIS");
%     scatter(Z_cart{nTrack}(1,:),Z_cart{nTrack}(3,:),75,'g',"o",'DisplayName',"Radar");
% end
% ylabel('$y$','Interpreter', 'Latex');
% xlabel('$x$','Interpreter', 'Latex');
% legend;
% hold off;
% allAxes = findall(figure(fc), 'Type', 'axes');
% for k = 1:length(allAxes)
%     allAxes(k).FontSize = 14;
% end

for nTrack = 1:N
    fc = fc + 1; figure(fc);
    %track = tracksTab(nTrack);
    hold on;
    scatter(Y_cart{nTrack}(:,1),Y_cart{nTrack}(:,2),75,'r',"+",'DisplayName',"AIS");
    scatter(Z_cart{nTrack}(:,1),Z_cart{nTrack}(:,2),75,'g',"o",'DisplayName',"Radar");
    ylabel('$y$','Interpreter', 'Latex');
    xlabel('$x$','Interpreter', 'Latex');
    title(nTrack);
    lgd = legend;
    lgd.FontSize = 16;
    hold off;
    allAxes = findall(figure(fc), 'Type', 'axes');
    for k = 1:length(allAxes)
        allAxes(k).FontSize = 14;
    end
end


for nTrack = 1:length(tracksTab)
    fc = fc + 1; figure(fc);
    %subplot(3,2,1);
    track = tracksTab(nTrack);
    %plot(estX{track}(:),'LineWidth',estValueLine.size);
    % plot(FiltX{track}(:,2),'LineWidth',estValueLine.size,'DisplayName','$\alpha_t(2)$');
    % hold on;
    % plot(filtXMean{track},'LineWidth',estValueLine.size,'DisplayName','Averaged $\alpha_t(2)$')
    % plot(xMarkov{track}(:),'LineWidth',refValueLine.size,'Color',refValueLine.color,'DisplayName','$x_t$');
    % plot(0.75.*ones(1,T),'m','LineWidth',refValueLine.size,'DisplayName','Threshold');
    % xlabel('$t$', 'Interpreter', 'Latex');
    % legend('Interpreter','latex');
    % hold off;
    % subplot(3,2,2);
    % track = tracksTab(nTrack);
    % %plot(estX{track}(:),'LineWidth',estValueLine.size);
    % plot(FiltX{track}(:,2),'LineWidth',estValueLine.size,'DisplayName','$\alpha_t(2)$');
    % hold on;
    % plot(xMarkov{track}(:),'LineWidth',refValueLine.size,'Color',refValueLine.color,'DisplayName','$x_t$');
    % xlabel('$t$', 'Interpreter', 'Latex');
    % legend('Interpreter','latex');
    % hold off;
    subplot(2,2,1);
    plot(Y_cart{track}(:,1),'LineWidth',estValueLine.size,'DisplayName',"AIS");
    hold on;
    plot(Z_cart{track}(:,1),'LineWidth',estValueLine.size,'DisplayName',"Radar");
    ylabel('$x$','Interpreter', 'Latex');
    xlabel('$t$', 'Interpreter', 'Latex');
    legend;
    hold off;
    subplot(3,2,2);
    plot(Y{track}(:,2),'LineWidth',estValueLine.size,'DisplayName',"AIS");
    hold on;
    plot(Z{track}(:,2),'LineWidth',estValueLine.size,'DisplayName',"Radar");
    ylabel('$\theta$','Interpreter', 'Latex');
    xlabel('$t$', 'Interpreter', 'Latex');
    legend;
    hold off;
    subplot(2,2,3);
    plot(Y_cart{track}(:,2),'LineWidth',estValueLine.size,'DisplayName',"AIS");
    hold on;
    plot(Z_cart{track}(:,2),'LineWidth',estValueLine.size,'DisplayName',"Radar");
    ylabel('$y$','Interpreter', 'Latex');
    xlabel('$t$', 'Interpreter', 'Latex');
    legend;
    hold off;
    subplot(2,2,4);
    plot(Y{track}(:,1),'LineWidth',estValueLine.size,'DisplayName',"AIS");
    hold on;
    plot(Z{track}(:,1),'LineWidth',estValueLine.size,'DisplayName',"Radar");
    ylabel('$r$','Interpreter', 'Latex');
    xlabel('$t$', 'Interpreter', 'Latex');
    legend;
    hold off;

    allAxes = findall(figure(fc), 'Type', 'axes');
    for k = 1:length(allAxes)
        allAxes(k).FontSize = 14;
    end
    %figTab(fc) = gcf;
end

fc = fc+1; figure(fc);
tracksToDisplay = [52,18,9,15,7,3];
%tracksToDisplay = [3,12,18];
tracksToDisplay = [49,23,40];
for nTrack = 1:length(tracksToDisplay)
    %subplot(3,2,1);
    track = tracksToDisplay(nTrack);
    %plot(estX{track}(:),'LineWidth',estValueLine.size);
    % plot(FiltX{track}(:,2),'LineWidth',estValueLine.size,'DisplayName','$\alpha_t(2)$');
    % hold on;
    % plot(filtXMean{track},'LineWidth',estValueLine.size,'DisplayName','Averaged $\alpha_t(2)$')
    % plot(xMarkov{track}(:),'LineWidth',refValueLine.size,'Color',refValueLine.color,'DisplayName','$x_t$');
    % plot(0.75.*ones(1,T),'m','LineWidth',refValueLine.size,'DisplayName','Threshold');
    % xlabel('$t$', 'Interpreter', 'Latex');
    % legend('Interpreter','latex');
    % hold off;
    % subplot(3,2,2);
    % track = tracksTab(nTrack);
    % %plot(estX{track}(:),'LineWidth',estValueLine.size);
    % plot(FiltX{track}(:,2),'LineWidth',estValueLine.size,'DisplayName','$\alpha_t(2)$');
    % hold on;
    % plot(xMarkov{track}(:),'LineWidth',refValueLine.size,'Color',refValueLine.color,'DisplayName','$x_t$');
    % xlabel('$t$', 'Interpreter', 'Latex');
    % legend('Interpreter','latex');
    % hold off;
    subplot(3,1,nTrack);
    scatter(Y_cart{track}(:,1), Y_cart{track}(:,2), 75, 'r', "+",'DisplayName',"AIS");
    hold on;
    scatter(Z_cart{track}(:,1), Z_cart{track}(:,2), 75, 'g', "o",'DisplayName',"Radar");
    hold off;
    ylabel('$y$','Interpreter', 'Latex');
    xlabel('$x$', 'Interpreter', 'Latex');
    %title(track);
    legend;
    hold off;

    allAxes = findall(figure(fc), 'Type', 'axes');
    for k = 1:length(allAxes)
        allAxes(k).FontSize = 15;
    end
    %figTab(fc) = gcf;
end

