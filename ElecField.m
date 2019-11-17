function ElecField(T,type,k)
%T is the succeeding time,type is the type of the motion of the charged
%particle,k is the ratio of the largest velocity to the speed of light(set
%to be 1).type should be chosen from 'harmonic'(harmonic oscillation),'line' and 'circle'.
xmin=-8;
ymin=-8;
xmax=8;
ymax=8;
dx=0.1;
dy=0.1;
dt=0.1;
Tpre=max((xmax-xmin),(ymax-ymin));%time needed to obtain a stable field before the motion begins
ntheta=36;%numbers of electric field lines
syms t;
switch type
    case 'harmonic'
    fx=cos(k*t);
    fy=0;
    case 'line'
    fx=k*t;
    fy=0;
    case 'circle'
    fx=cos(k*t);
    fy=sin(k*t);
end
xgrid=xmin:dx:xmax;
ygrid=ymin:dy:ymax;
nx=length(xgrid);
ny=length(ygrid);
auxt0=0:dt:T;
nt=length(auxt0);
tpre=0:dt:Tpre;
epsilon=tpre(end)-Tpre+dt;
tpre=tpre-epsilon;%tpre(end) might not be Tpre,so give a shift to tpre so that tpre(end)=Tpre-dt and t0(N)(will appear afterwards)=Tpre-dt and t0(N+1)=Tpre
N=length(tpre);
t0=[tpre,(auxt0+Tpre)];

vx0_symb=diff(fx,t);
vy0_symb=diff(fy,t);
ax0_symb=diff(vx0_symb,t);
ay0_symb=diff(vy0_symb,t);

%velocity and acceleration after the motion begins
auxx0=NewSubs(fx,t,auxt0);
auxy0=NewSubs(fy,t,auxt0);
auxvx0=NewSubs(vx0_symb,t,auxt0);
auxvy0=NewSubs(vy0_symb,t,auxt0);
auxax0=NewSubs(ax0_symb,t,auxt0);
auxay0=NewSubs(ay0_symb,t,auxt0);

%velocity and acceleration are all zero before the motion begins
%position before the motion begins is always the position when the motion
%is to start
vx0=[zeros(1,N),auxvx0];
vy0=[zeros(1,N),auxvy0];
ax0=[zeros(1,N),auxax0];
ay0=[zeros(1,N),auxay0];
x0=[auxx0(1)*ones(1,N),auxx0];
y0=[auxy0(1)*ones(1,N),auxy0];
%for start points of electric lines
shtheta=2*pi*(1:ntheta)/ntheta;%show_theta
sx0=5*dx*cos(shtheta);
sy0=5*dy*sin(shtheta);
% used if the edge of the electric particle is needed
% balltheta=2*pi*(1:96)/96;
% ballx=5*dx*cos(balltheta);
% bally=5*dy*sin(balltheta);

for it=1:nt
    %the delayed time must before current time,and the time whose influence
    %field has already gone out of the viewed region(determined by xgrid
    %and ygrid) is omitted to raise efficiency.Note the time interval 
    %corresbonding to N elements is the upper limit of the time needed for
    %a wavelet travelling all cross the viewed region.The speed of light is
    %set to be 1.
    xtemp=x0(it:N+it);
    ytemp=y0(it:N+it);
    ttemp=t0(it:N+it);
    %convert to 2-D array
    X=ones(ny,1)*xgrid;
    Y=ygrid'*ones(1,nx);
    %convert to 3-D array
    midX=reshape(X,nx*ny,1)*ones(1,N+1);
    midY=reshape(Y,nx*ny,1)*ones(1,N+1);
    newX=reshape(midX,ny,nx,N+1);
    newY=reshape(midY,ny,nx,N+1);
    clear midX;
    clear midY;
    newxtemp=reshape(ones(nx*ny,1)*xtemp,ny,nx,N+1);
    newytemp=reshape(ones(nx*ny,1)*ytemp,ny,nx,N+1);
    newttemp=reshape(ones(nx*ny,1)*ttemp,ny,nx,N+1);
    newt0=t0(it+N)*ones(ny,nx,N+1);
    %to determine the delayed time.All use 3-D array calculation instead of
    %element calculation to raise efficiency
    [~,idx]=min(abs( sqrt((newX-newxtemp).^2+(newY-newytemp).^2) - (newt0-newttemp) ),[],3);
    Rxstar=X-xtemp(idx);
    Rystar=Y-ytemp(idx);
    Rstar=sqrt(Rxstar.^2+Rystar.^2);
    idx0=(it-1)*ones(ny,nx)+idx;
    %calculate the electric field
    %first calculate parameters at delayed time
    vxstar=vx0(idx0);
    vystar=vy0(idx0);
    axstar=ax0(idx0);
    aystar=ay0(idx0);
    Sstar=Rstar-( Rxstar.*vxstar+Rystar.*vystar );
    auxvecx=Rxstar-Rstar.*vxstar;
    auxvecy=Rystar-Rstar.*vystar;
    prod1=1-vxstar.^2-vystar.^2;
    prod2=Rxstar.*axstar+Rystar.*aystar;
    prod3=Rxstar.*auxvecx+Rystar.*auxvecy;
    Sstar3=Sstar.^3;
    %the field
    Ex=( auxvecx.*(prod1+prod2) - prod3.*axstar )./Sstar3;
    Ey=( auxvecy.*(prod1+prod2) - prod3.*aystar )./Sstar3; 
    %calculate the position at current time
    x0it=x0(N+it);
    y0it=y0(N+it);
    %calculate the start points for streamline()
    sx=x0it+sx0;
    sy=y0it+sy0;
%     used if the edge of the electric particle is needed
%     plotballx=ballx+x0it;
%     plotbally=bally+y0it;
    
    streamline(xgrid,ygrid,Ex,Ey,sx,sy);
    axis equal
    axis([xmin xmax ymin ymax]);
    hold on
    plot(x0it,y0it,'g.','MarkerSize',50);
    % used if the edge of the electric particle is needed
    % plot(plotballx,plotbally,'r-');
    M(it)=getframe;
    cla;
end
movie(M,1,1/dt)
end