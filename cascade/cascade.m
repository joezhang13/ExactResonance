clear;
load resR2500r60Out.mat;

opt = 1;      %option: 1 for direct cascade; 2 for inverse cascade
uplow = 1;    %option: 1 for lower bound; 2 for upper bound
nitr = 5;
na = 0;       %number of active modes (to determine the end-of-loop criterion)
Kitr = K0;
Kmag = ((Kitr.*Kx).^2+(Kitr.*Ky).^2).^0.5;
Kact = (Kmag > r);
theta = linspace(0,2*pi);
xd = r*cos(theta);
yd = r*sin(theta);
if opt == 2
    Kact = (Kmag < r & Kmag > 0);
    xD = R*cos(theta);
    yD = R*sin(theta);
end

for itr = 1 : nitr
    %% "Saturated" angle resonance
    sactx = reshape(Kact.*Kx,(2*R)^2,1);
    sacty = reshape(Kact.*Ky,(2*R)^2,1);
    izero = (sactx==0)&(sacty==0);
    sactx(izero) = [];
    sacty(izero) = [];
    sact = [sactx, sacty];

    if uplow == 2
        smag = (sactx.^2+sacty.^2).^0.5;
        smag = unique(smag);
        sact = [];
        for i = 1 : length(smag)
            s = smag(i);
            for x = 0 : ceil(s)
                for y = 0 : ceil(s)
                    if norm([x,y]) == s
                        temp = [x, y; x, -y; -x, y; -x, -y];
                        temp = unique(temp, 'rows');
                        sact = [sact; temp];
                        for j = 1 : size(temp,1)
                            ikx = temp(j,1) + R + 1;
                            iky = temp(j,2) + R + 1;
                            Kact(ikx, iky) = 1;
                            Kitr(ikx, iky) = 1;
                        end
                    end
                end
            end
        end

        figure;
        plot(sact(:,1), sact(:,2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
        hold on
        axis equal;
        xlim([-R R]);
        ylim([-R R]);
        xlabel('$k_x$','Interpreter','latex');
        ylabel('$k_y$','Interpreter','latex');
        set(gca,'FontSize',15);

        smag = (sact(:,1).^2+sact(:,2).^2).^0.5;
        smag = unique(smag);
        if opt == 1
            kmax = max(smag);
            xD = kmax*cos(theta);
            yD = kmax*sin(theta);
            plot(xd, yd, '--b', 'LineWidth', 1.2);
            plot(xD, yD, '--c', 'LineWidth', 1.2);
            fill(xd, yd, 'r');
        else
            kmin = min(smag);
            xk = kmin*cos(theta);
            yk = kmin*sin(theta);
            plot(xd, yd, '--b', 'LineWidth', 1.2);
            plot(xk, yk, '--c', 'LineWidth', 1.2);
            fill([xd,xD], [yd,yD], 'r');
            plot([xd(end) xD(1)], [yd(end) yD(1)], '-r');
        end
%         saveas(gcf, ['itr', num2str(itr), 'Angle.fig']);
    end
    
    %% Generate modes from scale resonances involving new modes
    if opt == 1
        nkmagD = zeros(size(kmagD,1), 4);    %number of modes in (0~r,r~kmax,kmax~R,R~)
        nkmagD(:,1) = sum(kmagD<=r, 2);
        nkmagD(:,2) = sum(kmagD>r & kmagD<=kmax, 2);
        nkmagD(:,3) = sum(kmagD>kmax & kmagD<=R, 2);
        nkmagD(:,4) = sum(kmagD>R, 2);
        idx = nkmagD(:,1)<3 & nkmagD(:,3)<2 & nkmagD(:,4)==0;
    else
        nkmagD = zeros(size(kmagD,1), 4);    %number of modes in (0~kmin,kmin~r,r~R,R~)
        nkmagD(:,1) = sum(kmagD<kmin, 2);
        nkmagD(:,2) = sum(kmagD>=kmin & kmagD<r, 2);
        nkmagD(:,3) = sum(kmagD>=r & kmagD<=R, 2);
        nkmagD(:,4) = sum(kmagD>R, 2);
        idx = nkmagD(:,1)<2 & nkmagD(:,3)<3 & nkmagD(:,4)==0;
    end
    Q = qtD(idx, :);
    nkmagD = nkmagD(idx, :);

    for i = 1 : size(Q, 1)
        q = Q(i, :);
        nq = nkmagD(i, :);
        if opt == 1              %criterion determing the number of new modes we need to consider for each quartet
            nqc = 3 - nq(1);    
        else
            nqc = 3 - nq(3);
        end
        nqi = 0;
        for j = 1 : 4
            ktemp = q(2*j-1 : 2*j);
            findk = ktemp(1)==sact(:,1) & ktemp(2)==sact(:,2);
            if sum(findk) == 1
                nqi = nqi + 1;
            end
        end
        if nqi == nqc
            for j = 1 : 4
                idx = q(2*j-1) + R + 1;
                idy = q(2*j) + R + 1;
                Kitr(idx, idy) = 1;
            end
        end
    end
    Kmag = ((Kitr.*Kx).^2+(Kitr.*Ky).^2).^0.5;
    if opt == 1
        Kact = (Kmag > r);
    else
        Kact = (Kmag < r & Kmag > 0);
    end
    
    sx = reshape(Kact.*Kx,(2*R)^2,1);
    sy = reshape(Kact.*Ky,(2*R)^2,1);
    izero = (sx==0)&(sy==0);
    sx(izero) = [];
    sy(izero) = [];
    figure;
    plot(sx, sy, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
    hold on
    axis equal;
    xlim([-R R]);
    ylim([-R R]);
    xlabel('$k_x$','Interpreter','latex');
    ylabel('$k_y$','Interpreter','latex');
    set(gca,'FontSize',15);
    
    smag = (sx.^2+sy.^2).^0.5;
    smag = unique(smag);
    if opt == 1
        kmax = max(smag);
        xD = kmax*cos(theta);
        yD = kmax*sin(theta);
        plot(xd, yd, '--b', 'LineWidth', 1.2);
        plot(xD, yD, '--c', 'LineWidth', 1.2);
        fill(xd, yd, 'r');
    else
        kmin = min(smag);
        xk = kmin*cos(theta);
        yk = kmin*sin(theta);
        plot(xd, yd, '--b', 'LineWidth', 1.2);
        plot(xk, yk, '--c', 'LineWidth', 1.2);
        fill([xd,xD], [yd,yD], 'r');
        plot([xd(end) xD(1)], [yd(end) yD(1)], '-r');
    end
    %Plot the distribution at each iteraction
%     if uplow == 2
%         saveas(gcf, ['itr', num2str(itr), 'Scale.fig']);
%     else
%         saveas(gcf, ['itr', num2str(itr), '.fig']);
%     end
    %End-of-loop criterion
    if na == length(sx)
        break;
    else
        na = length(sx);
    end
end

if opt == 1
    if uplow == 1
        save(['kmaxR',num2str(R),'r',num2str(r),'Lower.mat'], 'R', 'r', 'kmax');
    else
        save(['kmaxR',num2str(R),'r',num2str(r),'Upper.mat'], 'R', 'r', 'kmax');
    end
else
    if uplow == 1
        save(['kminR',num2str(R),'r',num2str(r),'Lower.mat'], 'R', 'r', 'kmin');
    else
        save(['kminR',num2str(R),'r',num2str(r),'Upper.mat'], 'R', 'r', 'kmin');
    end
end
