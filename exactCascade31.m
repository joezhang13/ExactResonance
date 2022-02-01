clear;
load resN50.mat;

N=450;
kx=-N/2:N/2-1;
ky=-N/2:N/2-1;
[Kx,Ky]=meshgrid(kx,ky);
Kmag = (Kx.^2+Ky.^2).^0.5;
% dk = 0.00001;    %Broadening length for k
% dw = abs(dk./(2*(Kx.^2+Ky.^2).^0.25));   %Broadening length for w
dw = 10^(-9); %Broadening length for w
f = @(k1,k2,k3) ((k1(1)+k2(1)-k3(1)).^2+(k1(2)+k2(2)-k3(2)).^2).^0.25+...
                 (k3(1).^2+k3(2).^2).^0.25-(k1(1)^2+k1(2)^2)^0.25-(k2(1)^2+k2(2)^2)^0.25;

% K0 = (Kmag>=1)&(Kmag<=40); %Initial modes
kx0 = reshape(K0.*Kx,N^2,1);
ky0 = reshape(K0.*Ky,N^2,1);
iz0 = (kx0==0)&(ky0==0);
kx0(iz0) = [];
ky0(iz0) = [];
% figure;
% subplot(2,2,1);
% plot(kx0, ky0, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
% axis equal;
% xlim([kx(1) kx(end)]);
% ylim([ky(1) ky(end)]);
% xlabel('$k_x$','Interpreter','latex');
% ylabel('$k_y$','Interpreter','latex');
% set(gca,'FontSize',15);

Nitr = 3;
solf0 = K0;
solf = zeros(N, N, Nitr+1);
solf(:,:,1) = solf0;
for itr = 1 : Nitr
    kxx = reshape(solf(:,:,itr).*Kx,N^2,1);
    kyy = reshape(solf(:,:,itr).*Ky,N^2,1);
    iz = (kxx==0)&(kyy==0);
    kxx(iz) = [];
    kyy(iz) = [];
    kxky = [kxx,kyy];
    m = length(kxx);
    
    solf(:,:,itr+1) = solf(:,:,itr);
    for i = 1 : m
        for j = i : m
            for k = j : m
                dwf123 = f(kxky(i,:),kxky(j,:),kxky(k,:));
                dwf231 = f(kxky(j,:),kxky(k,:),kxky(i,:));
                dwf312 = f(kxky(k,:),kxky(i,:),kxky(j,:));
                if abs(dwf123) < dw
                    kk = kxky(i,:)+kxky(j,:)-kxky(k,:);
                    if (max(abs(kk)) < N/2) %&& ~isequal(kk, kxky(k,:))
                        ikk = [kk(2)+N/2+1, kk(1)+N/2+1];
                        solf(ikk(1),ikk(2),itr+1) = 1;
                    end
                end
                if abs(dwf231) < dw
                    kk = kxky(j,:)+kxky(k,:)-kxky(i,:);
                    if (max(abs(kk)) < N/2) %&& ~isequal(kk, kxky(i,:))
                        ikk = [kk(2)+N/2+1, kk(1)+N/2+1];
                        solf(ikk(1),ikk(2),itr+1) = 1;
                    end
                end
                if abs(dwf312) < dw
                    kk = kxky(k,:)+kxky(i,:)-kxky(j,:);
                    if (max(abs(kk)) < N/2) %&& ~isequal(kk, kxky(j,:))
                        ikk = [kk(2)+N/2+1, kk(1)+N/2+1];
                        solf(ikk(1),ikk(2),itr+1) = 1;
                    end
                end
            end
        end
    end
    
%     sx = reshape(solf(:,:,itr+1).*Kx,N^2,1);
%     sy = reshape(solf(:,:,itr+1).*Ky,N^2,1);
%     izero = (sx==0)&(sy==0);
%     sx(izero) = [];
%     sy(izero) = [];
%     subplot(2,2,itr+1);
%     plot(sx, sy, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
%     axis equal;
%     xlim([kx(1) kx(end)]);
%     ylim([ky(1) ky(end)]);
%     xlabel('$k_x$','Interpreter','latex');
%     ylabel('$k_y$','Interpreter','latex');
%     set(gca,'FontSize',15);
end

save('N50.mat', 'solf0', 'solf','Kx','Ky','N','kx','ky');
% saveas(gcf,'k9_11.fig');