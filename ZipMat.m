function mat = ZipMat(vec)

% This function takes an input column vector and zips it into a SQUARE matrix
c = sqrt(length(vec));
mat = zeros(c,c);
for(i=1:c)
    for(j=1:c)
        mat(i,j) = vec(i+(c*(j-1)),1);
    end
end


end