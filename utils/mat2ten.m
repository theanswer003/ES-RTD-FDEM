function T = mat2ten(mat, sz, mode)
N = length(sz);
newsz = [sz(mode), sz(1:mode-1), sz(mode+1:N)];
T = reshape(mat, newsz);
[~, idx] = sort([mode, 1:mode-1, mode+1:N]);
T = permute(T, idx);
