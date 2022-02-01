clear;
load resD2500d140Out.mat;

Kitr = K0;
Kmag = ((Kitr.*Kx).^2+(Kitr.*Ky).^2).^0.5;
Kout = (Kmag > d);
theta = linspace(0,2*pi);
xd = d*cos(theta);
yd = d*sin(theta);
for itr = 1 : 12
    %% "Saturated" angle resonance
    soutx = reshape(Kout.*Kx,N^2,1);
    souty = reshape(Kout.*Ky,N^2,1);
    izero = (soutx==0)&(souty==0);
    soutx(izero) = [];
    souty(izero) = [];
    sout = [soutx, souty];
%     smag = (soutx.^2+souty.^2).^0.5;
%     smag = unique(smag);
%     sout = [];
%     for i = 1 : length(smag)
%         s = smag(i);
%         for x = 0 : ceil(s)
%             for y = 0 : ceil(s)
%                 if norm([x,y]) == s
%                     temp = [x, y; x, -y; -x, y; -x, -y];
%                     temp = unique(temp, 'rows');
%                     sout = [sout; temp];
%                     for j = 1 : size(temp,1)
%                         ikx = temp(j,1) + N/2 + 1;
%                         iky = temp(j,2) + N/2 + 1;
%                         Kout(ikx, iky) = 1;
%                         Kitr(ikx, iky) = 1;
%                     end
%                 end
%             end
%         end
%     end
% 
%     figure;
%     plot(sout(:,1), sout(:,2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
%     hold on
%     axis equal;
%     xlim([-D D]);
%     ylim([-D D]);
%     xlabel('$k_x$','Interpreter','latex');
%     ylabel('$k_y$','Interpreter','latex');
%     set(gca,'FontSize',15);
%     
%     smag = (sout(:,1).^2+sout(:,2).^2).^0.5;
%     smag = unique(smag);
%     kmax = max(smag);
%     xD = kmax*cos(theta);
%     yD = kmax*sin(theta);
%     plot(xd, yd, '--b', 'LineWidth', 1.2);
%     plot(xD, yD, '--c', 'LineWidth', 1.2);
%     fill(xd, yd, 'r');
%     saveas(gcf, ['itr', num2str(itr), 'Angle.fig']);
    
    %% Generate modes from scale resonances involving new modes
    nkmagD = zeros(size(kmagD,1), 4);    %number of modes in (0~d,d~kmax,kmax~D,D~)
    nkmagD(:,1) = sum(kmagD<=d, 2);
    nkmagD(:,2) = sum(kmagD>d & kmagD<=kmax, 2);
    nkmagD(:,3) = sum(kmagD>kmax & kmagD<=D, 2);
    nkmagD(:,4) = sum(kmagD>D, 2);
    idx = nkmagD(:,1)<3 & nkmagD(:,3)<2 & nkmagD(:,4)==0;
    Q = qtD(idx, :);
    nkmagD = nkmagD(idx, :);

    for i = 1 : size(Q, 1)
        q = Q(i, :);
        nq = nkmagD(i, :);
        nqc = 3 - nq(1);    %criterion determing the number of new modes we need to consider for each quartet
        nqi = 0;
        for j = 1 : 4
            ktemp = q(2*j-1 : 2*j);
            findk = ktemp(1)==sout(:,1) & ktemp(2)==sout(:,2);
            if sum(findk) == 1
                nqi = nqi + 1;
            end
        end
        if nqi == nqc
            for j = 1 : 4
                idx = q(2*j-1) + N/2 + 1;
                idy = q(2*j) + N/2 + 1;
                Kitr(idx, idy) = 1;
            end
        end
    end
    Kmag = ((Kitr.*Kx).^2+(Kitr.*Ky).^2).^0.5;
    Kout = (Kmag > d);
    
    sx = reshape(Kout.*Kx,N^2,1);
    sy = reshape(Kout.*Ky,N^2,1);
    izero = (sx==0)&(sy==0);
    sx(izero) = [];
    sy(izero) = [];
    figure;
    plot(sx, sy, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
    hold on
    axis equal;
    xlim([-D D]);
    ylim([-D D]);
    xlabel('$k_x$','Interpreter','latex');
    ylabel('$k_y$','Interpreter','latex');
    set(gca,'FontSize',15);
    
    smag = (sx.^2+sy.^2).^0.5;
    smag = unique(smag);
    kmax = max(smag);
    xD = kmax*cos(theta);
    yD = kmax*sin(theta);
    plot(xd, yd, '--b', 'LineWidth', 1.2);
    plot(xD, yD, '--c', 'LineWidth', 1.2);
    fill(xd, yd, 'r');
%     saveas(gcf, ['itr', num2str(itr), 'Scale.fig']);
end