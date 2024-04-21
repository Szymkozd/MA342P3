function vec = unZipMat(mat)

% This function takes an input matrix and unzips it into a column matrix
[r,c] = size(mat);
vec = zeros(r*c,1);
for(i=1:r)
    for(j=1:c)
        vec(i+(r*(j-1)),1)=mat(i,j);
    end
end

end