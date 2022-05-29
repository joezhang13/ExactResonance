clear;
load kmaxDataR2500.mat;

figure;
plot(r, kmaxUpper, '-ro', 'LineWidth', 1.2);
hold on
plot(r, kmaxLower, '-bs', 'LineWidth', 1.2);
xlabel('$r_D$','Interpreter','latex');
ylabel('$k_{max}$','Interpreter','latex');
%legend('Upper limit', 'Lower limit');
set(gca,'Fontsize',15);