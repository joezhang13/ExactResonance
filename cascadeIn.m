clear;
load resD1100d500In.mat;

Kitr = K0;
Kmag = ((Kitr.*Kx).^2+(Kitr.*Ky).^2).^0.5;
Kin = (Kmag < d & Kmag > 0);
theta = linspace(0,2*pi);
xd = d*cos(theta);
yd = d*sin(theta);
xD = D*cos(theta);
yD = D*sin(theta);
for itr = 1 : 8
    %% "Saturated" angle resonance
    sinx = reshape(Kin.*Kx,(2*D)^2,1);
    siny = reshape(Kin.*Ky,(2*D)^2,1);
    izero = (sinx==0)&(siny==0);
    sinx(izero) = [];
    siny(izero) = [];
    sIn = [sinx, siny];
%     smag = (sinx.^2+siny.^2).^0.5;
%     smag = unique(smag);
%     sIn = [];
%     for i = 1 : length(smag)
%         s = smag(i);
%         for x = 0 : ceil(s)
%             for y = 0 : ceil(s)
%                 if norm([x,y]) == s
%                     temp = [x, y; x, -y; -x, y; -x, -y];
%                     temp = unique(temp, 'rows');
%                     sIn = [sIn; temp];
%                     for j = 1 : size(temp,1)
%                         ikx = temp(j,1) + D + 1;
%                         iky = temp(j,2) + D + 1;
%                         Kin(ikx, iky) = 1;
%                         Kitr(ikx, iky) = 1;
%                     end
%                 end
%             end
%         end
%     end
% 
%     figure;
%     plot(sIn(:,1), sIn(:,2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
%     hold on
%     axis equal;
%     xlim([-D D]);
%     ylim([-D D]);
%     xlabel('$k_x$','Interpreter','latex');
%     ylabel('$k_y$','Interpreter','latex');
%     set(gca,'FontSize',15);
%     
%     smag = (sIn(:,1).^2+sIn(:,2).^2).^0.5;
%     smag = unique(smag);
%     kmin = min(smag);
%     xk = kmin*cos(theta);
%     yk = kmin*sin(theta);
%     plot(xd, yd, '--b', 'LineWidth', 1.2);
%     plot(xk, yk, '--c', 'LineWidth', 1.2);
%     fill([xd,xD], [yd,yD], 'r');
%     plot([xd(end) xD(1)], [yd(end) yD(1)], '-r');
%     saveas(gcf, ['itr', num2str(itr), 'Angle.fig']);
    
    %% Generate modes from scale resonances involving new modes
    nkmagD = zeros(size(kmagD,1), 4);    %number of modes in (0~kmin,kmin~d,d~D,D~)
    nkmagD(:,1) = sum(kmagD<kmin, 2);
    nkmagD(:,2) = sum(kmagD>=kmin & kmagD<d, 2);
    nkmagD(:,3) = sum(kmagD>=d & kmagD<=D, 2);
    nkmagD(:,4) = sum(kmagD>D, 2);
    idx = nkmagD(:,1)<2 & nkmagD(:,3)<3 & nkmagD(:,4)==0;
    Q = qtD(idx, :);
    nkmagD = nkmagD(idx, :);

    for i = 1 : size(Q, 1)
        q = Q(i, :);
        nq = nkmagD(i, :);
        nqc = 3 - nq(3);    %criterion determing the number of new modes we need to consider for each quartet
        nqi = 0;
        for j = 1 : 4
            ktemp = q(2*j-1 : 2*j);
            findk = ktemp(1)==sIn(:,1) & ktemp(2)==sIn(:,2);
            if sum(findk) == 1
                nqi = nqi + 1;
            end
        end
        if nqi == nqc
            for j = 1 : 4
                idx = q(2*j-1) + D + 1;
                idy = q(2*j) + D + 1;
                Kitr(idx, idy) = 1;
            end
        end
    end
    Kmag = ((Kitr.*Kx).^2+(Kitr.*Ky).^2).^0.5;
    Kin = (Kmag < d & Kmag > 0);
    
    sinx = reshape(Kin.*Kx,(2*D)^2,1);
    siny = reshape(Kin.*Ky,(2*D)^2,1);
    izero = (sinx==0)&(siny==0);
    sinx(izero) = [];
    siny(izero) = [];
    figure;
    plot(sinx, siny, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
    hold on
    axis equal;
    xlim([-D D]);
    ylim([-D D]);
    xlabel('$k_x$','Interpreter','latex');
    ylabel('$k_y$','Interpreter','latex');
    set(gca,'FontSize',15);
    
    smag = (sinx.^2+siny.^2).^0.5;
    smag = unique(smag);
    kmin = min(smag);
    xk = kmin*cos(theta);
    yk = kmin*sin(theta);
    plot(xd, yd, '--b', 'LineWidth', 1.2);
    plot(xk, yk, '--c', 'LineWidth', 1.2);
    fill([xd,xD], [yd,yD], 'r');
    plot([xd(end) xD(1)], [yd(end) yD(1)], '-r');
%     saveas(gcf, ['itr', num2str(itr), 'Scale.fig']);
end