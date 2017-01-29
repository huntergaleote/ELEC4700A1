%% ELEC 4700 Assignment 1
%% Question 1
%Various constants are defined
%Boltzman's Constant
kb=1.38064852*10^(-23);
%Average Temperature
T=300;
%Effective electron mass
me=0.26*9.10938356*10^(-31);
%Here the thermal velocity is calculated
vth=sqrt(2*kb*T/(me));
maxx=200*10^(-9);
maxy=100*10^(-9);
tmn=0.2*10^(-12);
%Calculation for mean free path
MFP=vth*tmn;
e=ones(1000,4);
Px=ones(1000,1000);
Py=ones(1000,1000);
%Here uniform random values are assigned to x, y, and Vx
%Vy is calculated to ensure the total velocity is Vth
for i=1:1000
    x=rand*maxx;
    y=rand*maxy;
    vx=rand*vth*(-1)^floor(2*rand);
    vy=sqrt(vth^2-vx^2)*(-1)^floor(2*rand);
    e(i,1)=x;
    e(i,2)=y;
    e(i,3)=vx;
    e(i,4)=vy;
end
dt=100*10^(-9)/(vth*100);
%This is the main loop, each electron's position is changed 
%according to it's speed and then the boundary conditions 
%are checked and applied
figure(1)
title('Trajectory plot for Part 1');
for c=1:1000
    dx=dt.*e(:,3);
    dy=dt.*e(:,4);
    newx=e(:,1)+dx;
    newy=e(:,2)+dy;
    %Reversing Vy when the electron bounces off the top 
    %or bottom
    e(:,4)=e(:,4).*(-1).^((newy>maxy*ones(1000,1))+(newy<zeros(1000,1)));
    %Placing the electron at the opposite side if it exits
    %through the left or right
    e(:,1)=newx-maxx.*(newx>maxx*ones(1000,1))+maxx.*(newx<zeros(1000,1));
    %Reflecting electrons off the top and bottom
    e(:,2)=newy-2*(newy-maxy).*(newy>maxy*ones(1000,1))-2*newy.*(newy<zeros(1000,1));
    %Storing position data for plot
    Px(c,:)=e(:,1);
    Py(c,:)=e(:,2);
    plot(Px(1:c,1),Py(1:c,1),'.');
    %Plotting
    hold on
    for i=1:20
        plot(Px(1:c,i),Py(1:c,i),'.');
    end
    %Calculating and displaying current average temperature
    Tav(c)=mean(sqrt(e(:,3).^2+e(:,4).^2))^2*me/(2*kb);
    strtemp=['Average Temp: ',num2str(Tav(c))];
    text(0,0.5*10^(-8),strtemp)
    hold off
end
%Plotting average temperature
figure(2)
plot(Tav)
title('Temperature plot for Part 1');
xlabel('Cycle number');
ylabel('Average Temperature(K)');
%% Question 2
%Defining mean velocity and assigning velocity values in a
%normal disttribution
mu=sqrt(vth^2/2);
for i=1:1000
    x=rand*maxx;
    y=rand*maxy;
    vx=mu*randn;
    vy=mu*randn;
    e(i,1)=x;
    e(i,2)=y;
    e(i,3)=vx;
    e(i,4)=vy;
end
Pscat=zeros(1000,1);
DTscat=zeros(1000,1);
figure(3)
title('Trajectory plot for Part 2');
hold off
for c=1:1000
    dx=dt.*e(:,3);
    dy=dt.*e(:,4);
    newx=e(:,1)+dx;
    newy=e(:,2)+dy;
    %Calculating scattering probability and generating 
    %random values to see if electron scatters. If the
    %electron scatters the time since last scatter is reset
    Pscat=1-exp(-DTscat./tmn);
    Scat=rand(1000,1);
    DTscat=(DTscat+dt).*(Pscat<Scat);
    %Resetting Vx and Vy when scattering occurs
    e(:,3)=e(:,3).*(Pscat<Scat)+(randn.*mu.*ones(1000,1)).*(Pscat>Scat);
    e(:,4)=e(:,4).*(-1).^((newy>maxy*ones(1000,1))+(newy<zeros(1000,1)));
    e(:,4)=e(:,4).*(Pscat<Scat)+(randn.*mu.*ones(1000,1)).*(Pscat>Scat);
    e(:,1)=newx-maxx.*(newx>maxx*ones(1000,1))+maxx.*(newx<zeros(1000,1));
    e(:,2)=newy-2*(newy-maxy).*(newy>maxy*ones(1000,1))-2*newy.*(newy<zeros(1000,1));
    Px(c,:)=e(:,1);
    Py(c,:)=e(:,2);
    plot(Px(1:c,1),Py(1:c,1),'.');
    hold on
    for i=1:20
        plot(Px(1:c,i),Py(1:c,i),'.');
    end
    Tav(c)=mean(sqrt(e(:,3).^2+e(:,4).^2).^2.*me./(2.*kb));
    strtemp=['Average Temp: ',num2str(Tav(c))];
    text(0,0.5*10^(-8),strtemp)
    hold off
end
figure(4)
plot(Tav)
title('Temperature plot for Part 2');
xlabel('Cycle number');
ylabel('Average Temperature(K)');
%Creating histogram for final temperatures
figure(5)
Tf=sqrt(e(:,3).^2+e(:,4).^2).^2.*me./(2.*kb);
histogram(Tf);
title('Temperature Histogram for Part 2');
%% Question 3
%Defining the boundaries of the boxes
bx1=ones(1000,1)*0.8*10^(-7);
bx2=ones(1000,1)*1.2*10^(-7);
by1=ones(1000,1)*0.4*10^(-7);
by2=ones(1000,1)*0.6*10^(-7);
boxx=[bx1(1),bx1(1),bx2(1),bx2(1)];
box1y=[0,by1(1),by1(1),0];
box2y=[maxy,by2(1),by2(1),maxy];
%This can be set to 1 for specular boundaries or zero 
%for diffusive 
specular=1;
%Random variable setting loops untill electron is outside
%of the boxes
for i=1:1000
    inside=true;
    while inside;
        x=rand*maxx;
        y=rand*maxy;
        vx=mu*randn;
        vy=mu*randn;
        inside=(x>bx1)&(x<bx2)&((y<by1)|(y>by2));
        e(i,1)=x;
        e(i,2)=y;
        e(i,3)=vx;
        e(i,4)=vy;
    end
end
Pscat=zeros(1000,1);
DTscat=zeros(1000,1);
figure(6)
title('Trajectory plot for Part 3');
hold off
for c=1:1000
    dx=dt.*e(:,3);
    dy=dt.*e(:,4);
    newx=e(:,1)+dx;
    newy=e(:,2)+dy;
    Pscat=1-exp(-DTscat./tmn);
    Scat=rand(1000,1);
    DTscat=(DTscat+dt).*(Pscat<Scat);
    %Determining whether the electrons are bouncing off in 
    %the x or y directions
    xbound=((newx>bx1).*(e(:,1)<bx1)+(newx<bx2).*(e(:,1)>bx2)).*((newy<by1)+(newy>by2));
    ybound=((newy<by1).*(e(:,2)>by1)+(newy>by2).*(e(:,2)<by2)).*((newx>bx1).*(newx<bx2));
    %Extra terms added for handeling collisions with the boxes
    %for both types of barriers
    e(:,3)=e(:,3).*(Pscat<Scat)+(randn.*mu.*ones(1000,1)).*(Pscat>Scat);
    e(:,3)=e(:,3).*(-1).^(xbound*specular);
    e(:,3)=e(:,3).*(specular|(~xbound)|(~ybound))+abs(randn.*mu.*ones(1000,1)).*((~specular).*(xbound|ybound)).*(-1).^(((e(:,1)<bx1)&xbound)|(floor(2*rand).*ybound));
    e(:,4)=e(:,4).*(-1).^(ybound*specular);
    e(:,4)=e(:,4).*(specular|(~ybound)|(~xbound))+abs(randn.*mu.*ones(1000,1)).*((~specular).*(xbound|ybound)).*(-1).^(((newy>by2)&ybound)|(floor(2*rand).*xbound));
    e(:,4)=e(:,4).*(-1).^((newy>maxy*ones(1000,1))+(newy<zeros(1000,1)));
    e(:,4)=e(:,4).*(Pscat<Scat)+(randn.*mu.*ones(1000,1)).*(Pscat>Scat);
    newx=newx.*(~xbound)+(xbound).*((e(:,1)<bx1).*(2.*bx1-newx)+(e(:,1)>bx2).*(2.*bx2-newx));
    newy=newy.*(~ybound)+(ybound).*((newy<by1).*(2.*by1-newy)+(newy>by2).*(2.*by2-newy));
    e(:,1)=newx-maxx.*(newx>maxx*ones(1000,1))+maxx.*(newx<zeros(1000,1));
    e(:,2)=newy-2*(newy-maxy).*(newy>maxy*ones(1000,1))-2*newy.*(newy<zeros(1000,1));
    Px(c,:)=e(:,1);
    Py(c,:)=e(:,2);
    plot(Px(1:c,1),Py(1:c,1),'.');
    hold on
    for i=1:20
        plot(Px(1:c,i),Py(1:c,i),'.');
    end
    plot(boxx,box1y,'k');
    plot(boxx,box2y,'k');
    Tav(c)=mean(sqrt(e(:,3).^2+e(:,4).^2).^2.*me./(2.*kb));
    strtemp=['Average Temp: ',num2str(Tav(c))];
    text(0,0.5*10^(-8),strtemp)
    hold off
end
Tf=sqrt(e(:,3).^2+e(:,4).^2).^2.*me./(2.*kb);
%Creating an electron density map
figure(7)
for i=1:1000
 fx(i)=e(i,1);
 fy(i)=e(i,2);
end
title('Electron Density Map');
dmap=hist3([fy',fx'],[20,20]);
pcolor(dmap);


