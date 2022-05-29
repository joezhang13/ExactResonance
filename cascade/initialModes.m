clear;
load resR2500.mat;

R = 2500;
r = 60;
opt = 1;    %options: 1 for direct cascade, 2 for inverse cascade

kx = -R:R-1;
ky = -R:R-1;
[Kx,Ky] = meshgrid(kx,ky);
[K0, qtd, qtD, kmagd, kmagD, kmax, kmin] = initialize(r, R, sol, ncls, ng, opt);

%Plot the initial modes and the first iteration
sx = reshape(K0.*Kx,(2*R)^2,1);
sy = reshape(K0.*Ky,(2*R)^2,1);
izero = (sx==0)&(sy==0);
sx(izero) = [];
sy(izero) = [];
figure;
plot(sx, sy, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
hold on;
axis equal;
if opt == 1
    xlim([-2*kmax 2*kmax]);
    ylim([-2*kmax 2*kmax]);
else
    xlim([-R R]);
    ylim([-R R]);
end
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
set(gca,'FontSize',15);

theta = linspace(0,2*pi);
xd = r*cos(theta);
yd = r*sin(theta);
if opt == 1
    xk = kmax*cos(theta);
    yk = kmax*sin(theta);
else
    xk = kmin*cos(theta);
    yk = kmin*sin(theta);
end
xD = R*cos(theta);
yD = R*sin(theta);
plot(xd, yd, '-k', 'LineWidth', 1.2);
plot(xk, yk, '--c', 'LineWidth', 1.2);
if opt == 1
    fill(xd, yd, 'r');    %direct
else
    fill([xd,xD], [yd,yD], 'r');    %inverse
end
% plot(xD, yD, '-k', 'LineWidth', 1.2);
xlim([-R R]);
ylim([-R R]);   
xD = R*cos(theta);
yD = R*sin(theta);
plot(xd, yd, '-k', 'LineWidth', 1.2);
plot(xk, yk, '--c', 'LineWidth', 1.2);
if opt == 1
    fill(xd, yd, 'r');    %direct
else
    fill([xd,xD], [yd,yD], 'r');    %inverse
end
% plot(xD, yD, '-k', 'LineWidth', 1.2);
xlim([-R R]);
ylim([-R R]);

if opt == 1
    save(['resR',num2str(R),'r',num2str(r),'Out.mat'],... 
        'K0', 'qtd', 'qtD', 'kmagd', 'kmagD', 'Kx', 'Ky', 'kx', 'ky', 'R', 'r', 'kmax');
else
    save(['resR',num2str(R),'r',num2str(r),'In.mat'],... 
        'K0', 'qtd', 'qtD', 'kmagd', 'kmagD', 'Kx', 'Ky', 'kx', 'ky', 'R', 'r', 'kmin');
end
save(['quartR',num2str(R),'.mat'], 'R', 'Kx', 'Ky', 'kx', 'ky', 'kmagD', 'qtD');