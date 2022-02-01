clear;
load resN5000.mat;

D = 2500;
d = 200;

kx = -D:D-1;
ky = -D:D-1;
[Kx,Ky] = meshgrid(kx,ky);
opt = 1;    %options: 1 for direct cascade, 2 for inverse cascade
[K0, qtd, qtD, kmagd, kmagD, kmax, kmin] = initialize(d, D, sol, ncls, ng, opt);

%Plot the initial modes and the first iteration
sx = reshape(K0.*Kx,(2*D)^2,1);
sy = reshape(K0.*Ky,(2*D)^2,1);
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
    xlim([-D D]);
    ylim([-D D]);
end
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
set(gca,'FontSize',15);

theta = linspace(0,2*pi);
xd = d*cos(theta);
yd = d*sin(theta);
if opt == 1
    xk = kmax*cos(theta);
    yk = kmax*sin(theta);
else
    xk = kmin*cos(theta);
    yk = kmin*sin(theta);
end
xD = D*cos(theta);
yD = D*sin(theta);
plot(xd, yd, '--b', 'LineWidth', 1.2);
plot(xk, yk, '--c', 'LineWidth', 1.2);
plot(xD, yD, '-k', 'LineWidth', 1.2);

save('resD2500d200Out.mat', 'K0', 'qtd', 'qtD', 'kmagd', 'kmagD', 'Kx', 'Ky', 'kx', 'ky', 'N', 'D', 'd', 'kmax');
% save('resD1100d500In.mat', 'K0', 'qtd', 'qtD', 'kmagd', 'kmagD', 'Kx', 'Ky', 'kx', 'ky', 'N', 'D', 'd', 'kmin');