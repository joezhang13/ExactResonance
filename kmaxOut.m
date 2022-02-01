clear;
load kmaxOutDataD2500.mat;

figure;
plot(d, kmaxUpper, '-ro', 'LineWidth', 1.2);
hold on
plot(d, kmaxLower, '-bs', 'LineWidth', 1.2);
xlabel('$r$','Interpreter','latex');
ylabel('$k_{max}$','Interpreter','latex');
legend('Upper limit', 'Lower limit');
set(gca,'Fontsize',15);