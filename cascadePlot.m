clear;
load k200.mat;

n = size(solf, 3);
for i = 1 : n
    sx = reshape(solf(:,:,i).*Kx,N^2,1);
    sy = reshape(solf(:,:,i).*Ky,N^2,1);
    izero = (sx==0)&(sy==0);
    sx(izero) = [];
    sy(izero) = [];
    figure;
    plot(sx, sy, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
    axis equal;
    xlim([-N/2 N/2-1]);
    ylim([-N/2 N/2-1]);
    xlabel('$k_x$','Interpreter','latex');
    ylabel('$k_y$','Interpreter','latex');
    set(gca,'FontSize',15);
end