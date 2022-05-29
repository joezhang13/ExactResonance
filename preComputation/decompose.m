function [ArDg, nDec] = decompose(m)

sqm = m/floor(m^0.5);
sq = (1:sqm).^2;
mm = m./sq;
idxmm = (mm==fix(mm));
mm = mm(idxmm);
sq = sq(idxmm);
lm = length(mm);

ArDg = [];
for i = 1 : lm    %this loop is set in case that m is a multiple of a integer square and the sum of 2 squares (m=a^2(mm^2+nn^2))
    mi = mm(i);
    sqi = sq(i);
    mt = ceil((mi^2/4+1)/mi) - 1;    %maximum possible number of t
    t = ((1:mt)*mi-1).^0.5;
    t = t((t==fix(t))&(t>0));        %if t is a positive integer, keep it
    nt = length(t);
    
    if nt ~= 0
        xytemp = zeros(2, nt*2);
        for j = 1 : nt
            r0 = mi;
            r1 = t(j);
            r2 = r0 - floor(r0/r1)*r1;
            while r1^2 >= mi
                r0 = r1;
                r1 = r2;
                r2 = r0 - floor(r0/r1)*r1;
            end
            xytemp(:, j) = [r1; r2];        %[x;y] (x>y)
            xytemp(:, j+nt) = [r2; r1];     %[x;y] (x<y)
        end
        ArDg = [ArDg, xytemp*sqi^0.5];
    end
end
                
x2 = m^0.5;
xyeq = (m/2)^0.5;
if x2 == fix(x2)
    xytemp = [x2, 0; 0, x2];
    ArDg = [ArDg, xytemp];
end
if xyeq==fix(xyeq)
    xytemp = [xyeq; xyeq];
    ArDg = [ArDg, xytemp];
end

nDec = size(ArDg, 2);    %maximal number of decompositions

end