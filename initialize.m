function [K0, qtd, qtD, kmagd, kmagD, kmax, kmin] = initialize(d, D, sol, ncls, ng, opt)

%Choose the quartets with at least 3 wave vectors in the initial region (|k|<d or d<|k|<D)
K0 = zeros(2*D, 2*D);
qtd = [];    %quartets with at least 3 modes in the initial region
qtD = [];    %quartets in the whole domain (|k|<D)
for i = 1 : ncls
    for j = 1 : ng(i)
        if size(sol{i}{j},1) ~= 0
            ksol = sol{i}{j};
            kxx = ksol(:, [1,3,5,7]);
            kyy = ksol(:, [2,4,6,8]);
            kmag = (kxx.^2+kyy.^2).^0.5;
            nsD = sum(kmag<D, 2);    %number of wave vectors in |k|<D
            tempD = ksol(nsD==4, :);    %quartets with all wave vectors in |k|<D
            qtD = [qtD; tempD];
%             qtD = [qtD; ksol];
            if opt == 1
                nsd = sum(kmag<d, 2);    %number of wave vectors in |k|<d
            else
%                 nsd = sum(kmag>d & kmag<D, 2);    %number of wave vectors in d<|k|<D
                nsd = sum(kmag>d, 2);    %number of wave vectors in d<|k|<D
            end
            tempd = ksol(nsd>=3 & nsD==4, :);    %each quartet should have at least 3 wave vectors in |k|<d or d<|k|<D
%             tempd = ksol(nsd>=3, :);    %each quartet should have at least 3 wave vectors in |k|<d or d<|k|<D
            nsol = size(tempd, 1);
            qtd = [qtd; tempd];
            nik = nsol*4;
            kxx = reshape(tempd(:, [1,3,5,7]), nik, 1);
            kyy = reshape(tempd(:, [2,4,6,8]), nik, 1);
            kmag = (kxx.^2+kyy.^2).^0.5;
            ikx = kxx(kmag < D)+D+1;
            iky = kyy(kmag < D)+D+1;
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