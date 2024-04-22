clear all
close all
clc

BC = 60;
interior_cond = 60;
exteriot_cond = 0;

constants = [429,235,10500];%k,c,rho

constants(1) = constants(1)/100;%convert from meters to cm
constants(3) = constants(3)/1000000;

n = 100;%dimension of all matricies
nworker = 4;%set amount of workers for parallel for loop

[x,y,ztype,zrate] = slotFun(10,10,n);%set initial map 

deltax = x(2,1)-x(1,1);%use deltax and the calculated vaues to set deltat to a stable value using the equation in the textbook

deltat = constants(2)* constants(3)*(deltax^4)/(2*constants(1)*(2*(deltax^2)));
tf = 2;

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

job_size = floor(length(POI_INT)/nworker);%size of parallel job to be set to each individual worker

query_POI = zeros([nworker,job_size+mod(length(POI_INT),job_size)]);

for(i=1:(nworker-1))
    query_POI(i,1:job_size) = POI_INT(1,((i-1)*job_size)+1:i*job_size)-1;%segment points of interest into groups for each individual worker
end
query_POI(nworker,:)=POI_INT(1,((nworker-1)*job_size)+1:end)-1;%shove the remaining points into last worker


for(i=1:steps)
    zout = zeros([nworker,length(z)]);%reset zout
    
    parfor(quad=1:nworker)
        modifier = 0;
        if(quad==nworker)
            modifier = length(query_POI)-job_size;
        end
    zout(quad,:) = PDEsolvePARALLEL(x,y,z,ztype,zrate,deltat,constants,query_POI(quad,1:(job_size+modifier)));%fill zout in parallel
    end
    z = sum(zout,1);%combine resulting zout values to form z
    z(POI_BC_DIRC)=BC;%set exterior points and boundaries to the specified conditions

    z(z<0)=0;%set temps equal to 0 when they get low enough
    surf(ZipMat(x')-scaling,ZipMat(y')-scaling,ZipMat(z')-scaling)
    drawnow
    M(i) = getframe;%store step as a frame

    if(any(isnan(z))||any(isinf(z)))
        fprintf("%f|%f\n",sum(isnan(z)),sum(isinf(z)))
        find(isnan(z));
        error("OUCH")
        break%if a value comes up as nan or inf, halt
    end
end

fig.Visible = 'on';

movie(M,1);