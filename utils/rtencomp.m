function [B, Qlist] = rtencomp(X, R, p, q)
N = length(R);
B = X;
% p = 10;
% q = 2;
Qlist = cell(N-1, 1);


for mode = 1:N
    Bn = ten2mat(B, mode);
    [m, n] = size(Bn);
    l = R(mode) + p;
    k = min(m, l);
    a = -1.0; b = 1.0;
    O = a + (b - a).*rand(n, k);
    Y = Bn * O;
    if q > 1
        for i = 1:q
            [Y, ~] = lu(Y);
            [Z, ~] = lu(Bn' * Y);
            Y = Bn * Z;
        end
    end
    [Qn, ~] = qr(Y, 0);
    B = ttm(B, Qn', mode);
    Qlist{mode} = Qn;
end