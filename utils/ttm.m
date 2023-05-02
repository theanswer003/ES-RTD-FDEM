function newT = ttm(T, mat, mode)
sz = size(T);
newsz = sz;
newsz(mode) = size(mat, 1);
Tn = ten2mat(T, mode);
newT = mat*Tn;
newT = mat2ten(newT, newsz, mode);