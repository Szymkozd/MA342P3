function [x,y,type,rate] = slotFun(D,L,n)

% slotFun constructs an nxn map of a slot shapped object with a width of D
% and a side length (only the straight portion) of L
% Outputs : 
% x,y = coordinates of point
% type = what this point is (Neummann BC:3, Diriclet BC:2, interior:1, exterior:0)

x = zeros(n,n);%pre allocated memory for positions and types
y = zeros(n,n);
type = zeros(n,n);
rate = zeros(n,n);

nbc = -1;

step = (L+D)/(n-2);%set the step size so that the box will be slightly larger than the slot is tall

scale = 1/2;%scaling factor to assist with the thickness of the boundaries

tol = scale/step;%ceil(log10(n^(2/3)));%determines generally how thick the boundaries will be
county = 1;
countx = 1;
%nrange = (-(L+D)/2-step*tol:step:(L+D)/2)+step*tol;
nrange = linspace(-(L+D)/2-step*tol,(L+D)/2+step*tol,n);
for(i=1:n)
    
    for(j=1:n)
        x(countx,county) = nrange(j);
        y(countx,county) = nrange(i);
        if(abs(x(countx,county))>D/2+tol*step)%set all values to the left and right of the slot to outside points
            current_type = 0;
        elseif(abs(y(countx,county))<L/2&&abs(x(countx,county))<D/2-tol*step)%sets all values within center "square section" to inside points
            current_type=1;
        elseif((D/2)^2>((x(countx,county))^2)+((abs(y(countx,county))-L/2)^2)&&abs(y(countx,county))>L/2)%sets all values within hemispheres to inside points
            current_type = 1;
        elseif((D/2+tol*step)^2<((abs(x(countx,county)))^2)+((abs(y(countx,county))-L/2)^2)&&abs(y(countx,county))>L/2)%sets all values outside hemisphere to outside points
            current_type = 0;
        elseif(abs(y(countx,county))<L/2&&abs(x(countx,county))<D/2)%set square bc's
            current_type = 2;
        else
            current_type = 3;%the only remaining points are the semi-circle boundaries
        end
        type(countx,county) = current_type;
        countx = countx + 1;
    end
    county = county + 1;
    countx = 1;
end

rate(type==3)=nbc;%set rate conditions

end