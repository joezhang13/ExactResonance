function [K0, qtd, qtD, kmagd, kmagD, kmax, kmin] = initialize(r, R, sol, ncls, ng, opt)

%Choose the quartets with at least 3 wave vectors in the initial region (|k|<r or r<|k|<R)
K0 = zeros(2*R, 2*R);
qtd = [];    %quartets with at least 3 modes in the initial region
qtD = [];    %quartets in the whole domain (|k|<R)
for i = 1 : ncls
    for j = 1 : ng(i)
        if size(sol{i}{j},1) ~= 0
            ksol = sol{i}{j};
            kxx = ksol(:, [1,3,5,7]);
            kyy = ksol(:, [2,4,6,8]);
            kmag = (kxx.^2+kyy.^2).^0.5;
            nsD = sum(kmag<R, 2);    %number of wave vectors in |k|<R
            tempD = ksol(nsD==4, :);    %quartets with all wave vectors in |k|<R
            qtD = [qtD; tempD];
%             qtD = [qtD; ksol];
            if opt == 1
                nsd = sum(kmag<r, 2);    %number of wave vectors in |k|<r
            else
%                 nsd = sum(kmag>d & kmag<D, 2);    %number of wave vectors in d<|k|<R
                nsd = sum(kmag>r, 2);    %number of wave vectors in r<|k|<R
            end
            tempd = ksol(nsd>=3 & nsD==4, :);    %each quartet should have at least 3 wave vectors in |k|<r or r<|k|<R
%             tempd = ksol(nsd>=3, :);    %each quartet should have at least 3 wave vectors in |k|<r or r<|k|<R
            nsol = size(tempd, 1);
            qtd = [qtd; tempd];
            nik = nsol*4;
            kxx = reshape(tempd(:, [1,3,5,7]), nik, 1);
            kyy = reshape(tempd(:, [2,4,6,8]), nik, 1);
            kmag = (kxx.^2+kyy.^2).^0.5;
            ikx = kxx(kmag < R)+R+1;
            iky = kyy(kmag < R)+R+1;
            for k = 1 : length(ikx)
                K0(ikx(k), iky(k)) = 1;
            end
        end
    end
end

%Calculate the magnitudes of modes in each quartet
kmagd = zeros(size(qtd,1),4);
kmagD = zeros(size(qtD,1),4);
for i = 1 : 4
    kmagd(:, i) = (qtd(:,2*i-1).^2 + qtd(:,2*i).^2).^0.5;
    kmagD(:, i) = (qtD(:,2*i-1).^2 + qtD(:,2*i).^2).^0.5;
end
if opt == 1
    kmax = max(max(kmagd));
    kmin = [];
else
    kmin = min(min(kmagd));
    kmax = [];
end

end