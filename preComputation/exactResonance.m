clear;
R = 2500;    %maximum wavenumber in the domain
 
%% Step 1: Class indexes
Arq = ones(1, floor(R^2));
pls = primes((R^2)^0.25);
%First pass
for p =  pls
    kappa = floor(R^2/p^4);
    idx = (1:kappa)*p^4;
    Arq(idx) = 0;
end
%Second pass
pls = primes(R^2);
puls = pls(mod(pls,4)==3);
for p = puls
    amax = floor(R^2/p);
    a = 1:amax;
    idx = a(mod(a,p)~=0)*p;
    Arq(idx) = 0;
    idx = a((mod(a,p)==0)&(mod(a,p^2)==0))*p;
    Arq(idx) = 0;
end
Wq = 1 : floor(R^2);
Wq = Wq(Arq==1);
Mq = floor((R^2./Wq).^0.25);    %maximum possible weight for a given q
piD = length(Wq);               %length of the class list

%% Step 2: Solving the sum-of-weight equation
keep = (Mq > 2);       %discard the classes with Mq = 1,2 which contain no solution
Wq = Wq(keep);
Mq = Mq(keep);
ncls = length(Wq);      %number of valid classes

sow = cell(1, ncls);    %solution of weights
for i = 1 : ncls
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
    lab = ones(ns, 1);    %label of cases (1,2,3,4)
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

%% Step 3: Decomposition into sum of two squares
ArD = cell(1, ncls);      %sets of decompositions
ArDtoW = cell(1, ncls);   %number of decompositions for each gamma^4*q
for i = 1 : ncls
    q = Wq(i);
    ArD{i} = cell(Mq(i), 1);
    Decq = zeros(Mq(i), 1);
    for gamma = 1 : Mq(i)
        m = gamma^4*q;
        [ArDg, nDec] = decompose(m);
        ArD{i}{gamma} = ArDg;
        Decq(gamma) = nDec;
    end
    ArDtoW{i} = Decq;
end

%% Step 4: Checking linear conditions
signs = ff2n(8);
signs = 2*signs-1;
sol = cell(1, ncls);    %set of solutions
ng = zeros(1, ncls);    %number of solutions
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

save(['resR', num2str(R), '.mat'], 'R', 'sol', 'ncls', 'ng');
