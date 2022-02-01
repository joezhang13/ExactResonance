clear;

N = 500;
D = N/2;
 
%% Step 1: Class indexes
Arq = ones(1, floor(2*D^2));
pls = primes((2*D^2)^0.25);
%First pass
for p =  pls
    kappa = floor(2*D^2/p^4);
    idx = (1:kappa)*p^4;
    Arq(idx) = 0;
end
%Second pass
pls = primes(2*D^2);
puls = pls(mod(pls,4)==3);
for p = puls
    amax = floor(2*D^2/p);
    a = 1:amax;
    idx = a(mod(a,p)~=0)*p;
    Arq(idx) = 0;
    idx = a((mod(a,p)==0)&(mod(a,p^2)==0))*p;
    Arq(idx) = 0;
%     if p <= (2*D^2)^0.5
%         idx = a(mod(a,p)~=0)*p;
%         Arq(idx) = 0;
%         if p <= (2*D^2)^(1/3)
%             idx = a((mod(a,p)==0)&(mod(a,p^2)==0))*p;
%             Arq(idx) = 0;
%         end
%     end
end
%Third pass
Wq = 1 : floor(2*D^2);
Wq = Wq(Arq==1);
Mq = floor((2*D^2./Wq).^0.25);
piD = length(Wq);

%% Step 2: Decomposition into sum of two squares
ArD = cell(1, piD);
ArDtoW = cell(1, piD);
for i = 1 : piD
    q = Wq(i);
    ArD{i} = cell(Mq(i), 1);
    Decq = zeros(Mq(i), 1);
    for gamma = 1 : Mq(i)
        m = gamma^4*q;
        [ArDg, nDec] = decompose(m, D);
%         mt = ceil((m^2/4+1)/m) - 1;    %maximum possible number of t
%         t = ((1:mt)*m-1).^0.5;
%         t = t((t==fix(t))&(t>0));           %if t is a positive integer, keep it
%         nt = length(t);
%         
% %         if nt == 0                        %this corresponds to the case where 
% %             sqm = m/floor(m^0.5);
% %             sq = (1:sqm).^2;
% %             mm = m./sq;
% %             idxmm = (mm==fix(mm));
% %             mm = mm(idxmm);
% %             sq = sq(idxmm);
% %             m = mm(end);
% %             sqr = sq(end);
% %             mt = ceil((m^2/4+1)/m) - 1;    %maximum possible number of t
% %             t = ((1:mt)*m-1).^0.5;
% %             t = t((t==fix(t))&(t>0));           %if t is a positive integer, keep it
% %             nt = length(t);
% %         end
%         
%         x2 = m^0.5;
%         xyeq = (m/2)^0.5;
%         if x2==fix(x2)
%             ArDg = zeros(2, nt*2+2);
%             ArDg(:, end-1) = [x2; 0];    %[x^2;0]
%             ArDg(:, end) = [0; x2];    %[0;x^2]
%         elseif xyeq==fix(xyeq)
%             ArDg = zeros(2, nt*2+1);
%             ArDg(:, end) = [xyeq; xyeq];    %[x;y] (x=y)
%         else
%             ArDg = zeros(2, nt*2);
%         end
%         Decq(gamma) = size(ArDg, 2);    %maximal number of decompositions
%         
%         for j = 1 : nt
%             r0 = m;
%             r1 = t(j);
%             r2 = r0 - floor(r0/r1)*r1;
%             while r1^2 >= m
%                 r0 = r1;
%                 r1 = r2;
%                 r2 = r0 - floor(r0/r1)*r1;
%             end
%             ArDg(:, j) = [r1; r2];        %[x;y] (x>y)
%             ArDg(:, j+nt) = [r2; r1];     %[x;y] (x<y)
%         end
        ArD{i}{gamma} = ArDg;
        Decq(gamma) = nDec;
%         if nDec == 0
%             ArD{i}(gamma) = [];
%             Decq(gamma) = [];
%             Mq(i) = Mq(i) - 1;
%         else
%             ArD{i}{gamma} = ArDg;
%             Decq(gamma) = nDec;
%         end
    end
    ArDtoW{i} = Decq;
end

%% Step 3: Solving the sum-of-weight equation
sow = cell(1, piD);
for i = 1 : piD
    ns = 0;
    gs = zeros(1, 4);
    for S = 2 : 2*Mq(i)
        if S <= Mq(i)+1
            for g1L = 1 : floor(S/2)
                g2L = S - g1L;
                for g1R = g1L : floor(S/2)
                    g2R = S - g1R;
                    ns = ns + 1;
                    gs(ns, :) = [g1L, g2L, g1R, g2R];
                end
            end
        else
            for g1L = (S-Mq(i)) : floor(S/2)
                g2L = S - g1L;
                for g1R = g1L : floor(S/2)
                    g2R = S - g1R;
                    ns = ns + 1;
                    gs(ns, :) = [g1L, g2L, g1R, g2R];
                end
            end
        end
    end
    lab = ones(ns, 1);
    for j = 1 : ns
        if gs(j, 1)==gs(j, 3)
            if gs(j, 3)==gs(j, 4)
                lab(j) = 4;
            else
                lab(j) = 2;
            end
        elseif gs(j, 3)==gs(j, 4)
            lab(j) = 3;
        end
    end
    gs = [gs, lab];
    sow{i} = gs;
end

%% Step 4: Discarding "lean" cases
lidx = ones(1, piD);
for i = 1 : piD
    if (Mq(i)==1) && (ArDtoW{i}<=4)
        lidx(i) = 0;
    end
end
Wq = Wq(lidx>0);
Mq = Mq(lidx>0);
ArD = ArD(lidx>0);
ArDtoW = ArDtoW(lidx>0);
sow = sow(lidx>0);
ncls = length(Wq);

%% Step 5: Checking linear conditions
signs = ff2n(8);
signs = 2*signs-1;
sol = cell(1, ncls);
ng = zeros(1, ncls);
for i = 1 : ncls
    gsol = sow{i};
    ng(i) = size(gsol, 1);
    sol{i} = cell(1, ng(i));
    for j = 1 : ng(i)
        g = gsol(j, 1:4);        % (g1L,g2L,g1R,g2R)
        cs = gsol(j, 5);         % case number
        
        % for scale resonances, we only need to consider case 1 and 3
        if (cs == 1) || (cs == 3)
            MN = zeros(1, 8);
            nDecq = ArDtoW{i}(g);    % number of decompositions
            mn = zeros(1, 8);        % (m1,n1,m2,n2,m3,n3,m4,n4)
            for j1 = 1 : nDecq(1)
                mn(1:2) = ArD{i}{g(1)}(:, j1)';
                for j2 = 1 : nDecq(2)
                    mn(3:4) = ArD{i}{g(2)}(:, j2)';
                    for j3 = 1 : nDecq(3)
                        mn(5:6) = ArD{i}{g(3)}(:, j3)';
                        for j4 = 1 : nDecq(4)
                            mn(7:8) = ArD{i}{g(4)}(:, j4)';
                            mnsol = mn.*signs;
                            idx = (mnsol(:,1)+mnsol(:,3))==(mnsol(:,5)+mnsol(:,7)) & ... 
                                  (mnsol(:,2)+mnsol(:,4))==(mnsol(:,6)+mnsol(:,8));
                            mnsol = mnsol(idx, :);
                            mnsol = unique(mnsol, 'rows');
                            if cs == 3
                                nls = size(mnsol, 1);
                                rmlb = zeros(nls, 1);
                                for k = 1 : nls
                                    temp = mnsol(k, 5:8);
                                    temp = temp([3,4,1,2]);
                                    for l = k + 1 : nls
                                        temp1 = mnsol(l, 5:8);
                                        if (rmlb(l)==0) && isequal(temp, temp1)
                                            rmlb(l) = 1;
                                        end
                                    end
                                end
                                mnsol = mnsol(rmlb==0, :);
                            end
                            MN = [MN; mnsol];
%                         is1 = (cs==1);
%                         is2 = (cs==2) & (mn(1)~=mn(5)) & (mn(3)~=mn(7));
%                         is3 = (cs==3) & (mn(5)~=mn(7));
%                         is4 = (cs==4) & (mn(1)~=mn(5)) & (mn(5)~=mn(7)) & (mn(7)~=mn(3));
%                         iscount = (is1 | is2 | is3 | is4);
%                         if iscount
%                             smL = mn(1) + mn(3);
%                             dmL = abs(mn(1) - mn(3));
%                             smR = mn(5) + mn(7);
%                             dmR = abs(mn(5) - mn(7));
%                         end
                        end
                    end
                end
            end
            if size(MN, 1) > 1
                MN = MN(2:end, :);
                sol{i}{j} = MN;
            end
        end
    end
end

% save('resN2000.mat', 'N', 'sol', 'ncls', 'ng');

%% Collecting the modes involved in exact resonances
% D = 250;
% d = 150;
% 
% kx = -D:D-1;
% ky = -D:D-1;
% [Kx,Ky] = meshgrid(kx,ky);
% opt = 2;    %options: 1 for direct cascade, 2 for inverse cascade
% [K0, qtd, qtD, kmagd, kmagD, kmax, kmin] = initialize(d, D, sol, ncls, ng, opt);
% 
% %Plot the initial modes and the first iteration
% sx = reshape(K0.*Kx,(2*D)^2,1);
% sy = reshape(K0.*Ky,(2*D)^2,1);
% izero = (sx==0)&(sy==0);
% sx(izero) = [];
% sy(izero) = [];
% figure;
% plot(sx, sy, 'ro', 'MarkerFaceColor','r', 'MarkerSize',2);
% hold on;
% axis equal;
% xlim([-D D]);
% ylim([-D D]);
% xlabel('$k_x$','Interpreter','latex');
% ylabel('$k_y$','Interpreter','latex');
% set(gca,'FontSize',15);
% 
% theta = linspace(0,2*pi);
% xd = d*cos(theta);
% yd = d*sin(theta);
% if opt == 1
%     xk = kmax*cos(theta);
%     yk = kmax*sin(theta);
% else
%     xk = kmin*cos(theta);
%     yk = kmin*sin(theta);
% end
% xD = D*cos(theta);
% yD = D*sin(theta);
% plot(xd, yd, '--b', 'LineWidth', 1.2);
% plot(xk, yk, '--c', 'LineWidth', 1.2);
% plot(xD, yD, '-k', 'LineWidth', 1.2);
% 
% % save('resD500d25Out.mat', 'K0', 'qtd', 'qtD', 'kmagd', 'kmagD', 'Kx', 'Ky', 'kx', 'ky', 'N', 'D', 'd', 'kmax');
% % save('resD200d150In.mat', 'K0', 'qtd', 'qtD', 'kmagd', 'kmagD', 'Kx', 'Ky', 'kx', 'ky', 'N', 'D', 'd', 'kmin');