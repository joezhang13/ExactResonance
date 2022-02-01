clear;
load kminInd250Data.mat;

figure;
plot(D, d-kminUpper, '-ro', 'LineWidth', 1.2);
hold on
plot(D, d-kminLower, '-bs', 'LineWidth', 1.2);
xlabel('$R$','Interpreter','latex');
ylabel('$r-k_{min}$','Interpreter','latex');
legend('Upper limit', 'Lower limit');
set(gca,'Fontsize',15);