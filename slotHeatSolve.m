clear all
close all
clc

BC = 60;
interior_cond = 60;
exteriot_cond = 0;

copper =    [ 401, 385,  8933];      %k,c,rho
nickel =    [90.7, 444,  8900];      %k,c,rho
tin =       [66.6, 227,  7310];      %k,c,rho
aluminum =  [ 237, 903,  2702];      %k,c,rho
silver =    [ 429, 235, 10500];      %k,c,rho
stainless = [15.1, 480,  8055];      %k,c,rho

constants = stainless;

constants(1) = constants(1)/100;%convert from meters to cm
constants(3) = constants(3)/1000000;

n = 100;%dimension of all matricies

[x,y,ztype,zrate] = slotFun(10,10,n);%set initial map 

deltax = x(2,1)-x(1,1);%use deltax and the calculated vaues to set deltat to a stable value using the equation in the textbook

deltat = constants(2)* constants(3)*(deltax^4)/(2*constants(1)*(2*(deltax^2)));
tf = 500;

steps = ceil(tf/deltat);%determine steps based on deltat and desired tf value, round to an integer

x = unZipMat(x)';%convert matricies to row vectors
y = unZipMat(y)';
ztype = unZipMat(ztype)';
z = ztype;%split into the solution matrix z and the type matrix to reference later on if needed
%find which indicies indicate boundary conditions, interior points, and exterior points
POI_BC_NEUM = find(z==3);%neumann boundary conditions
POI_BC_DIRC = find(z==2);%Dirichlet bondary conditions
POI_INT = find(z==1);%interior
POI_EXT = find(z==0);%exterior

%set values to initial conditions

z(POI_BC_DIRC) = BC;
z(POI_INT) = interior_cond;
z(POI_EXT) = exteriot_cond;
z(POI_BC_NEUM) = exteriot_cond;%neumann conditions are placed outside of the geometry

scaling = exteriot_cond;

fig = figure;
surf(ZipMat(x')-scaling,ZipMat(y')-scaling,ZipMat(z')-scaling)
%surf(x,y,z)
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';

M(steps) = struct('cdata',[],'colormap',[]);%initialize struct to hold movie frames

fig.Visible = 'off';%hide frames while calculating

for(i=1:steps)
    z = PDEsolve(x,y,z,ztype,zrate,deltat,constants);
    z(POI_BC_DIRC)=BC;%set exterior points and boundaries to the specified conditions
    z(POI_BC_NEUM)=0;
    z(POI_EXT)=exteriot_cond;
    %z(POI_BC_NEUM)=mean(z(POI_BC_NEUM));
    z(z<0)=0;%set temps equal to 0 when they get low enough
    surf(ZipMat(x')-scaling,ZipMat(y')-scaling,ZipMat(z')-scaling)
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