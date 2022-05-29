clear;
load quartR2500.mat;

%% Collect data
dk = 20;                 %delta k used to calculate the average
KD = zeros(2*R, 2*R);    %indicator of active mode 
Kq = zeros(2*R, 2*R);    %multiplicity of each mode
nq = size(qtD, 1);       %number of quartets
Oq = zeros(2*R, 2*R);    %sum of order index
for i = 1 : nq
    for j = 1 : 4
        idx = qtD(i,2*j-1) + R + 1;
        idy = qtD(i,2*j) + R + 1;
        KD(idx,idy) = 1;
        Kq(idx,idy) = Kq(idx,idy) + 1;
        switch j
            case 1
                Oq(idx,idy) = Oq(idx,idy) + 1;
            case 2
                Oq(idx,idy) = Oq(idx,idy) + 4;
            case 3
                Oq(idx,idy) = Oq(idx,idy) + 2;
            case 4
                Oq(idx,idy) = Oq(idx,idy) + 3;
        end
    end
end
qx = reshape(KD.*Kx,(2*R)^2,1);
qy = reshape(KD.*Ky,(2*R)^2,1);
iz = (qx==0)&(qy==0);
qx(iz) = [];
qy(iz) = [];
mq = reshape(Kq,(2*R)^2,1);        %multiplicity at different scales
mq(iz) = [];
Kmag = (Kx.^2 + Ky.^2).^0.5;
Kmag = reshape(Kmag,(2*R)^2,1);
Kmag(iz) = [];
% Kmag = Kmag(Kmag<=R & Kmag>0);
% mq = mq(Kmag<=R & Kmag>0);

%% Plot Ms and Os
dkk = 50;
nn = 2*R/dkk;
Kqq = mat2cell(Kq, ones(1,nn)*dkk, ones(1,nn)*dkk);
MM = zeros(nn, nn);
for i = 1 : nn
    for j = 1 : nn
        Mij=Kqq{i,j};
        MM(i,j)=mean(mean(Mij(Mij~=0)));
    end
end
Kxx = mat2cell(Kx, ones(1,nn)*dkk, ones(1,nn)*dkk);
Kxc = cellfun(@mean,cellfun(@mean,Kxx,'UniformOutput',false));
Kyy = mat2cell(Ky, ones(1,nn)*dkk, ones(1,nn)*dkk);
Kyc = cellfun(@mean,cellfun(@mean,Kyy,'UniformOutput',false));
figure;
contourf(Kxc,Kyc,MM,100,'LineColor','none');
axis equal;
colorbar;
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
set(gca,'FontSize',15);

NOq = Kq;
NOq(NOq == 0) = 1;
Oq = Oq./NOq;        %order index calculated by the average of the rankings
Oqq = mat2cell(Oq, ones(1,nn)*dkk, ones(1,nn)*dkk);
CC = zeros(nn, nn);
for i = 1 : nn
    for j = 1 : nn
        Cij=Oqq{i,j};
        CC(i,j)=mean(mean(Cij(Cij~=0)));
    end
end
figure;
contourf(Kxc,Kyc,CC,100,'LineColor','none');
axis equal;
colorbar;
caxis([1 4]);
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
set(gca,'FontSize',15);

%% Compute angle-averaged Ms and Os
kmag = reshape(kmagD, nq*4, 1);
[kmag, l] = sort(kmag);     %multiplicity
kkmag = unique(kmag);       %scales
lambda = repmat([1,4,2,3], nq, 1);
lambda = reshape(lambda, nq*4, 1);
lambda = lambda(l);         %order of magnitude in each quartet (from small to large: 1,2,3,4)
nk = length(kkmag);
Nl = zeros(nk, 4);          %Nl(i,j): number of modes at a given scale kkmag(i) with the same order j
[kkmag, ia, ~] = unique(kmag);
ia = [ia; 4*nq+1];
for i = 1 : nk
    temp = lambda(ia(i):ia(i+1)-1);
    for j = 1 : 4
        Nl(i,j) = sum(temp==j);
    end
end
npar = R/dk;
kkp = linspace(0, R, npar+1);
kkp = kkp(1:end-1) + dk/2;
Stemp = zeros(npar,1);
ltemp = zeros(npar,4);
for i = 1 : nk
    ii = ceil(kkmag(i)/dk);
    Stemp(ii) = Stemp(ii) + sum(Nl(i, :));
    ltemp(ii, :) = ltemp(ii, :) + Nl(i, :);
end
MsTemp = zeros(npar, 1);
ntemp = zeros(npar, 1);
nm = length(mq);
for i = 1 : nm
    ii = ceil(Kmag(i)/dk);
    MsTemp(ii) = MsTemp(ii) + mq(i);
    ntemp(ii) = ntemp(ii) + 1;
end
Ms = MsTemp./ntemp;
Oqk = sum(ltemp.*[1,2,3,4], 2)./Stemp;

figure;
yyaxis left
plot(kkp, Ms, 'LineWidth', 1.2);
xlabel('$k$','Interpreter','latex');
ylabel('$\tilde{M}_s$','Interpreter','latex');
yyaxis right
plot(kkp, Oqk, 'LineWidth', 1.2);
ylabel('$\tilde{O}_s$','Interpreter','latex');
set(gca,'FontSize',15);
