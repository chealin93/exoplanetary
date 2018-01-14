%normalized by 1AU
figure(1);clf;
figure(2);clf;
figure(3);clf;
figure(4);clf;
figure(5);clf;

format long;
clear all;

global eps;
global m_t;
global number_of_planet;
global tot_num;

%-------------------------input parameter----------------------------------
%input function
%when you put a value in each input function, you must enter exact value
%not type space!!!!!
number_of_planet = 0; %input('put the number of additional planets');
yesno = 2; %input('do you want to generate random number of moons? answer yes(1) or no(2)');

if yesno == 2
    input_num_m = 0; %input('put the additional moon to each planets');
end

separation1 = 0.22431;        
separation2 = 0.7;   
inclination1 = 0;       
inclination2 = 0;       
%-----------------------------set exosystem--------------------------------
%recognize the exoplanetary system and properties based on solar system
%proportional expression (standard 'the erath')

%-------------------------basic planet setting-----------------------------
%m_s = first star, m_p = second star

G = 1;
n_p = number_of_planet + 2;
m_s = [0.6897];
m_p = [0.20255];
m_t = {m_s m_p};         %normalized by earth mass

se_t = {separation1};
inclination_t = {inclination1};

pos_s = [0 0 0];
dy1 = cosd(inclination1)*se_t{1};
dz1 = sind(inclination1)*se_t{1};
pos_p = [0 dy1 dz1];
pos_t = {pos_s pos_p};

x=0;y=0;z=0;mass=0;
for i=1:length(pos_t)
    mass = mass + m_t{i};
    x = x + m_t{i}*pos_t{i}(1);
    y = y + m_t{i}*pos_t{i}(2);
    z = z + m_t{i}*pos_t{i}(3);
end

x = x/mass;
y = y/mass;
z = z/mass;

CM_pos = [x y z];
CM_mass = mass;
    
vel_sx = sqrt(G*m_t{2}^2/(dist(pos_p,pos_s))/CM_mass);
vel_px = sqrt(G*m_t{1}^2/(dist(pos_s,pos_p))/CM_mass);

vel_s = [vel_sx 0 0];
vel_p = [-vel_px 0 0];
vel_t = {vel_s vel_p};

%-------------------------N-body simulation--------------------------------
%with fixxing the system of one star, one planet, one moon
%----------------------------planets system--------------------------------
%mass range of planet(0.5<=p<=10)
%up to jupiter mass
for i = 3:n_p
    m_t{i} = [0.5+(100-0.5).*rand];
end

for i = 2:n_p-1
    se_t{i} = 0.3+(1.5-0.3).*rand;
end

for i = 2:n_p-1
    inclination_t{i} = 0;%rand*90;
end

for i = 2:n_p-1
    dy1 = cosd(inclination_t{i})*se_t{i};
    dz1 = sind(inclination_t{i})*se_t{i};
    temp_pos = [0 dy1 dz1];
    pos_t{i+1}=temp_pos;
end

for i = 3:n_p
    temp_vel = sqrt(G*m_t{1}/dist(pos_t{1},pos_t{i}));
    vel_t{i} = [temp_vel 0 0];
end


%--------------------------basic moon setting------------------------------
m_c = [0.00033];
m_t{length(m_t)+1} = m_c;
se_t{length(se_t)+1}=separation2;
inclination_t{length(inclination_t)+1}=inclination2;

dy2 = cosd(inclination_t{end})*se_t{end};
dz2 = sind(inclination_t{end})*se_t{end};
a = CM_pos(2)+separation2;
b = separation2;
pos_c = [0 a dz2];

pos_t{length(pos_t)+1} = pos_c;

vel_cx = sqrt(G*CM_mass/dist(CM_pos,pos_c));
vel_c = [vel_cx 0 0];
index_moon = length(vel_t)+1;
vel_t{length(vel_t)+1} = vel_c;

%-----------------------additional moon setting----------------------------
num_list = [];
for i = 2:n_p
    if yesno == 1
        temp_num_m = randi([0,5]);
    else
        temp_num_m = input_num_m;
    end
    
    %random number of moons
    %temp_num_m = 2;  %for a moment, origin is 1
    %temp_num_m = randi([0,5]);
    num_list = [num_list temp_num_m];
    
    for j = 1:temp_num_m
        k = length(m_t) + 1;
        m_t{k} = [rand*0.05];
    end
   
    for j = 1:temp_num_m
        k = length(se_t) + 1;
        se_t{k} = 0.001+(0.005-0.001)*rand;
        inclination_t{k} = rand*90;
        if inclination_t{i-1}+inclination_t{k}>90
            ddy = cosd(inclination_t{i-1}-inclination_t{k})*se_t{k};
            ddz = sind(inclination_t{i-1}-inclination_t{k})*se_t{k};
        else
            ddy = cosd(inclination_t{i-1}+inclination_t{k})*se_t{k};
            ddz = sind(inclination_t{i-1}+inclination_t{k})*se_t{k};
        end
    
        temp_pos = pos_t{i}+[0 ddy ddz];
        pos_t{k+1} = temp_pos;
        
        temp_vel = vel_t{i} + [sqrt(G*m_t{i}/dist(pos_t{i},temp_pos)) 0 0];
        vel_t{k+1} = temp_vel;
    end
end

%----------------finding period of system approximately--------------------
% 
semi_a = dist(pos_t{end},CM_pos);
period = sqrt(4*pi^2/(m_t{3}+CM_mass)*semi_a^3);
%---------------------------integration part-------------------------------
m=1;eps=0;t0 =0;tf=100000;%2000000;
period_list = [];p = period;
while tf > p
    period_list = [period_list p];
    p = p + period; 
end

period_point_number = length(period_list)/10

tot_num = length(pos_t)
n = tot_num;

y0 = [pos_t{1:n} vel_t{1:n}];length(y0);
figure(1);hold on;
for pn = 1:period_point_number %(length(period_list)-1)
    options = odeset('RelTol',1e-10); %,'AbsTol',1e-12);
    [T,Y] = ode45('several_planet_moon_rk',[t0,period_list(pn)],y0,options);

    k=1;j=3*n+1;
    pos={};vel={};
    for i = 1:n
        pos{i}=Y(:,k:k+2);
        k=k+3;
        vel{i}=Y(:,j:j+2);
        j=j+3;
    end
    plot(pos{3}(end,1),pos{3}(end,2),'.');
    xlabel('x');
    ylabel('y');
    title('poincare map');
    
    t0 = period_list(pn);
    tf = period_list(pn+1);
    
    y0 = Y(end,:);
end

k=1;j=3*n+1;
pos={};vel={};
for i = 1:n
    pos{i}=Y(:,k:k+2);
    k=k+3;
    vel{i}=Y(:,j:j+2);
    j=j+3;
end


%nC2
clear D;l=length(T);
D=[];D_t={};q=1;
for i = 1:n-1
    for j = i+1:n
        for k =1:length(pos{1})
            D(k) = sqrt(sum((pos{i}(k,:)-pos{j}(k,:)).^2));
        end
        
        %relation of each object for Potential energy
        %separation and mass
        
        D_t{q} = D;
        R_m_t{q} = m_t{i}*m_t{j}; 
        q = q+1;
    end
end

%for figure 2 and 3 (poincare map) of star
m_vel = [];

%-------------------------------2d plot------------------------------------
figure(2);
%sun's x-position as a function of time
subplot(3,3,1);
    plot(T,pos{1}(:,1)-pos{end}(:,1));
    title('x-pos of sun');

%position of sun as function of time
subplot(3,3,2);
    plot(T,vel{1}(:,1));
    title('P');
    
%momentum of sun as function of time
subplot(3,3,3);
    plot(T,m_t{1}*vel{1}(:,1));
    title('Q');
    
%distance between two body
subplot(3,3,4);
    plot(T,D_t{1});hold on;
    title('distance between S&E');

%velocity
subplot(3,3,5);
    plot(T,vel{1},'g');hold on;
    plot(T,vel{2},'k');hold on;
%     plot(T,vel{index_moon},'r');
    axis([0 max(T)/10 -1000 1000]);
    title('velocity of S$E&M');

clear E U;L=length(T);

for i = 1:L
    E_t = 0;
    %kinetic energy
    for j = 1:n
        temp_E = 1/2*m_t{j}.*sum(vel{j}(i,:).^2);
        E_t = E_t + temp_E;
    end
    %potential energy
    U_t = 0;
    for j = 1:length(D_t)
        temp_U = -R_m_t{j}./D_t{j}(i);
        U_t = U_t +temp_U;
    end
    E(i) = E_t;
    U(i) = U_t;
    sval(i) = U(i)/E(i);
end

dT = diff(T);
eta(m) = max(dT)/min(dT)
eta_t  = dT/min(dT);

%eta as function of time
subplot(3,3,6);
    plot(T(1:end-1),eta_t);hold on;
    title('Eta as function of t');

%energy and error
U=U';E=E';
H = E + U;
  
subplot(3,3,7);
    plot(T,E);hold on
    title('kinetic');
subplot(3,3,8);
    plot(T,U);hold on
    title('potential');
subplot(3,3,9);
    plot(T,1-H/mean(H));hold on;
    title('1-tot/mean(tot)');

figure(3);
plot(pos{1}(:,1),vel{1}(:,1));
xlabel('pos');ylabel('vel');
title('projection of trajectories');

figure(4);
plot(vel{1}(:,1),m_t{1}*vel{1}(:,1));
xlabel('x');ylabel('y');
title('P-Q diagram');
%-------------------------------3d plot------------------------------------
figure(5);
plot(pos{1}(:,1),pos{1}(:,2),'r');hold on;
plot(pos{2}(:,1),pos{2}(:,2),'b');hold on;grid on;
plot(pos{index_moon}(:,1),pos{index_moon}(:,2),'k');grid on;
xlabel('x');ylabel('y');zlabel('z');
title('3d plot');

figure(6);
plot(T,E+U);
title('total energy');
