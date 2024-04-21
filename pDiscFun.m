function [x,y,type] = pDiscFun(p,n)

%%This function generates a map based on conditions from exercise 6 for a
%%given p

x = zeros(n,n);%pre allocated memory for positions and types
y = zeros(n,n);
type = zeros(n,n);

nbc = -1;

step = 1/(n-2);%set the step size so that the box will be slightly larger than the slot is tall

scale = 1/10;%scaling factor to assist with the thickness of the boundaries

tol = scale/step;%ceil(log10(n^(2/3)));%determines generally how thick the boundaries will be
county = 1;
countx = 1;

nrange = linspace(-1-step*tol,1+step*tol,n);

for(i=1:n)
    for(j=1:n)
        x(countx,county) = nrange(j);
        y(countx,county) = nrange(i);
        if(1<((abs(x(countx,county)^p)+abs(y(countx,county)^p))^(1/p)))%if outside circle set to outside points
            current_type = 0;
        elseif((1-tol*step)>((abs(x(countx,county)^p)+abs(y(countx,county)^p))^(1/p)))%set inside points
            current_type = 1;
        else%set dirichlet conditions
            current_type = 2;
        end
        type(countx,county) = current_type;
        countx = countx + 1;
    end
    county = county + 1;
    countx = 1;
end

end