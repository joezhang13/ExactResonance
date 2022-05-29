clear;
load kminDatar1000.mat;

figure;
plot(R, kminUpper, '-ro', 'LineWidth', 1.2);
hold on
plot(R, kminLower, '-bs', 'LineWidth', 1.2);
xlabel('$R$','Interpreter','latex');
ylabel('$k_{min}$','Interpreter','latex');
% legend('Upper limit', 'Lower limit');
set(gca,'Fontsize',15);
