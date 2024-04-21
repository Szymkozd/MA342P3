clear all
close all
clc

BC = 0.1;

constants = [2, 1, 0];


p=100;
n = 100;%dimension of all matricies

[x,y,ztype] = pDiscFun(p,n);%set initial map

deltax = x(2,1)-x(1,1);

deltat = 0.001;% keep at or below 0.0001
tf = 0.75;

steps = round(tf/deltat);

x = unZipMat(x)';%convert matricies to row vectors
y = unZipMat(y)';
ztype = unZipMat(ztype)';
z = ztype;%split into the solution matrix z and the type matrix to reference later on if needed
%find which indicies indicate boundary conditions, interior points, and exterior points
POI_BC_DIRC = find(z==2);%Dirichlet bondary conditions
POI_INT = find(z==1);%interior
POI_EXT = find(z==0);%exterior

%set values to initial conditions

z(POI_BC_DIRC) = BC;
z(POI_INT) = ((abs(x(POI_INT).^p)+abs(y(POI_INT).^p)).^(1./p))./10;%set interior conditions based on problem statement
zhist = z;%wave eqn needs history input

fig = figure;
surf(ZipMat(x'),ZipMat(y'),ZipMat(z'))

axis ([-1,1,-1,1,0,1])
ax = gca;
ax.NextPlot = 'replaceChildren';

M(steps) = struct('cdata',[],'colormap',[]);%initialize struct to hold movie frames

fig.Visible = 'off';%hide frames while calculating

for(i=1:steps)
    zhisthold = z;%hold onto the current value to use as history for the step after this
    z = PDEsolve(x,y,z,ztype,zhist,deltat,constants);
    zhist = zhisthold;%set new historical value
    z(POI_BC_DIRC)=BC;%set exterior points and boundaries to the specified conditions
    z(POI_EXT)=0;

    surf(ZipMat(x'),ZipMat(y'),ZipMat(z'))
    drawnow
    M(i) = getframe;%store step as a frame
    [bigz,bigi]=max(z);
    [smallz,smalli]=min(z);
    if(any(isnan(z))||any(isinf(z)))
        fprintf("%f|%f\n",sum(isnan(z)),sum(isinf(z)))
        find(isnan(z));
        error("OUCH")
        break%if a value comes up as nan or inf, halt
    end
end

fig.Visible = 'on';

movie(M,1,500);