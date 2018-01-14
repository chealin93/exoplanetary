%normalized by 1AU
close all;
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

separation1 = 1; %0.22431;        
separation2 = 0.7048;   
inclination1 = 0;       
inclination2 = 0;       

%-----------------------------set exosystem--------------------------------
%recognize the exoplanetary system and properties based on solar system
%proportional expression (standard 'the erath')

%-------------------------basic planet setting-----------------------------
%m_s = first star, m_p = second star

G = 1;
n_p = number_of_planet + 2;
m_s = [0.5];
m_p = [0.5];
m_t = {m_s m_p};         %normalized by earth mass

se_t = {separation1};
inclination_t = {inclination1};

pos_s = [0 0 0];
dy1 = cosd(inclination1)*se_t{1};
dz1 = sind(inclination1)*se_t{1};
pos_p = [0 dy1 dz1];
pos_t = {pos_s pos_p};

vel_sx = sqrt(G*m_t{2}^2/(dist(pos_p,pos_s)*(m_t{1}+m_t{2})));
vel_px = sqrt(G*m_t{1}^2/(dist(pos_s,pos_p)*(m_t{1}+m_t{2})));

vel_s = [vel_sx 0 0];
vel_p = [-0.9*vel_px 0 0];
vel_t = {vel_s vel_p};


x=0;y=0;z=0;
vx=0;vy=0;vz=0;
mass=0;
for i=1:length(pos_t)
    mass = mass + m_t{i};
    x = x + m_t{i}*pos_t{i}(1);
    vx = vx + m_t{i}*vel_t{i}(1);
    y = y + m_t{i}*pos_t{i}(2);
    vy = vy + m_t{i}*vel_t{i}(2);
    z = z + m_t{i}*pos_t{i}(3);
    vz = vz + m_t{i}*vel_t{i}(3);
end

x = x/mass;
y = y/mass;
z = z/mass;

vx = vx/mass;
vy = vy/mass;
vz = vz/mass;

CM_pos = [x y z];
CM_mass = mass;    
CM_vel = [vx vy vz];

vel_s = vel_s - CM_vel;
vel_p = vel_p - CM_vel;

%--check---- 
m_t{1}*vel_s + m_t{2}*vel_p;
%-----------


lumi_s = 0.148; %solar luminosity
lumi_p = 0.0057;

total_l = lumi_s + lumi_p;
r_i = sqrt(total_l/1.1);
r_o = sqrt(total_l/0.53);

%if two star is main sequence

%-------------------------N-body simulation--------------------------------
%with fixxing the system of one star, one planet, one moon
%----------------------------planets system--------------------------------
%mass range of planet(0.5<=p<=10)
%up to jupiter mass

% for i = 3:n_p
%     m_t{i} = [0.5+(100-0.5).*rand];
% end
% 
% for i = 2:n_p-1
%     se_t{i} = 0.3+(1.5-0.3).*rand;
% end
% 
% for i = 2:n_p-1
%     inclination_t{i} = 0;%rand*90;
% end
% 
% for i = 2:n_p-1
%     dy1 = cosd(inclination_t{i})*se_t{i};
%     dz1 = sind(inclination_t{i})*se_t{i};
%     temp_pos = [0 dy1 dz1];
%     pos_t{i+1}=temp_pos;
% end
% 
% for i = 3:n_p
%     temp_vel = sqrt(G*m_t{1}/dist(pos_t{1},pos_t{i}));
%     vel_t{i} = [temp_vel 0 0];
% end


%--------------------------basic moon setting------------------------------
% m_c = [0.00033];
% 
% m_t{length(m_t)+1} = m_c;
% se_t{length(se_t)+1}=separation2;
% inclination_t{length(inclination_t)+1}=inclination2;
% 
% dy2 = cosd(inclination_t{end})*se_t{end};
% dz2 = sind(inclination_t{end})*se_t{end};
% a = CM_pos(2)+separation2;
% b = separation2;
% pos_c = [0 a dz2]
% 
% pos_t{length(pos_t)+1} = pos_c;
% 
% vel_t{2} = planet velocity
% vel_cx = sqrt(G*CM_mass/dist(CM_pos,pos_c));
% vel_c = [vel_cx 0 0];
% index_moon = length(vel_t)+1;
% vel_t{length(vel_t)+1} = vel_c;

%-----------------------additional moon setting----------------------------
% num_list = [];
% for i = 2:n_p
%     if yesno == 1
%         temp_num_m = randi([0,5]);
%     else
%         temp_num_m = input_num_m;
%     end
%     
%     %random number of moons
%     %temp_num_m = 2;  %for a moment, origin is 1
%     %temp_num_m = randi([0,5]);
%     num_list = [num_list temp_num_m];
%     
%     for j = 1:temp_num_m
%         k = length(m_t) + 1;
%         m_t{k} = [rand*0.05];
%     end
%    
%     for j = 1:temp_num_m
%         k = length(se_t) + 1;
%         se_t{k} = 0.001+(0.005-0.001)*rand;
%         inclination_t{k} = rand*90;
%         if inclination_t{i-1}+inclination_t{k}>90
%             ddy = cosd(inclination_t{i-1}-inclination_t{k})*se_t{k};
%             ddz = sind(inclination_t{i-1}-inclination_t{k})*se_t{k};
%         else
%             ddy = cosd(inclination_t{i-1}+inclination_t{k})*se_t{k};
%             ddz = sind(inclination_t{i-1}+inclination_t{k})*se_t{k};
%         end
%     
%         temp_pos = pos_t{i}+[0 ddy ddz];
%         pos_t{k+1} = temp_pos;
%         
%         temp_vel = vel_t{i} + [sqrt(G*m_t{i}/dist(pos_t{i},temp_pos)) 0 0];
%         vel_t{k+1} = temp_vel;
%     end
% end

%----------------finding period of system approximately--------------------
% 
semi_a = dist(CM_pos,pos_t{2});
period = sqrt(4*pi^2/(CM_mass)*semi_a^3);

%--------------------------determining spt---------------------------------
m=1; eps=0; t0 =0; tf=100; 

bi_period = sqrt(4*pi^2*dist(pos_s,pos_p)^3/(G*(m_t{1}+m_t{2})));
spt = bi_period
spt_list = [];
while tf > spt
    spt_list = [spt_list spt];
    spt = spt + bi_period;
end

%---------------------------integration part-------------------------------

period_list = [];p = period;
while tf > p
    period_list = [period_list p];
    p = p + period;
end

period_point_number = length(period_list)

tot_num = length(pos_t)
n = tot_num;

y0 = [pos_t{1:n} vel_t{1:n}];length(y0);
options = odeset('RelTol',1e-12); %,'AbsTol',1e-12);
[T,Y] = ode45('several_planet_moon_rk',[t0,tf],y0,options);


k=1;j=3*n+1;
pos={};vel={};

for i = 1:n
    pos{i}=Y(:,k:k+2);
    k=k+3;
    vel{i}=Y(:,j:j+2);
    j=j+3;
end


%relation each other from 1 to end
l=length(T);
D=[];D_t={};q=1;
for i = 1:n-1
    for j = i+1:n
        for k =1:l
            D(k) = sqrt(sum((pos{i}(k,:)-pos{j}(k,:)).^2));
        end
                
        D_t{q} = D;
        R_m_t{q} = m_t{i}*m_t{j}; 
        q = q+1;
    end
end


%-------------------------------2d plot------------------------------------
figure('Name','Analysis','NumberTitle','off');
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
    title('distance of primaries');

%velocity
subplot(3,3,5);
    plot(T,vel{1},'g');hold on;
    plot(T,vel{2},'k');hold on;
%     plot(T,vel{index_moon},'r');
    axis([0 max(T)/10 -1000 1000]);
    title('velocity of S&E&M');

L=length(T);
E = []; U = []; TE = []; Ep = []; Up = []; TE_p = [];
for i = 1:L  
    
    E_1 = 0;
    temp_kie1 = 1/2*m_t{1}*sum(vel{1}(i,:).^2);
    E_1 = E_1 + temp_kie1;
    
    U_1 = 0;
    for j = 1:length(D_t)
        temp_poe1 = -R_m_t{j}/D_t{j}(i);
        U_1 = U_1 + temp_poe1;
    end
    
    E_2 = 0;
    temp_kie2 = 1/2*m_t{2}*sum(vel{2}(i,:).^2);
    E_2 = E_2 + temp_kie2;
    
    %secondary star potential energy
    U_2 = 0;
    for j = 1:length(D_t)
        temp_poe2 = -R_m_t{j}/D_t{j}(i);
        U_2 = U_2 + temp_poe2;
    end
    
    %total system's kinetic energy
    E_t = 0;
    for j = 1:length(pos_t)
        temp_kie = 1/2*m_t{j}*sum(vel{j}(i,:).^2);
        E_t = E_t + temp_kie;
    end
    
    %tptal system's potential energy
    U_t = 0;
    for j = 1:length(D_t)
        temp_poe = -R_m_t{j}./D_t{j}(i);
        U_t = U_t +temp_poe;
    end
   
    E(i) = E_t;
    U(i) = U_t;
    
    E1(i) = E_1;
    U1(i) = U_1;
    
    E2(i) = E_2;
    U2(i) = U_2;
    
end

E=E';U=U';E1=E1';U1=U1';E2=E2';U2=U2';

%---------------------interpolating requried values------------------------
interp_T = T(1):0.01*bi_period:T(end);
interp_T=interp_T';

interp_pos{1} = interp1(T,pos{1},interp_T);
interp_pos{2} = interp1(T,pos{2},interp_T);
interp_vel{1} = interp1(T,vel{1},interp_T);
interp_vel{2} = interp1(T,vel{2},interp_T);
interp_E = interp1(T,E,interp_T);
interp_U = interp1(T,U,interp_T);

spline_pos{1} = interp1(T,pos{1},interp_T);
spline_pos{2} = interp1(T,pos{2},interp_T);
spline_vel{1} = interp1(T,vel{1},interp_T);
spline_vel{2} = interp1(T,vel{2},interp_T);
spline_E = interp1(T,E,interp_T);
spline_U = interp1(T,U,interp_T);


%----------------------------exact solution--------------------------------
%spt is specific time, when you want to give the specific time,
%you must call like T(spt), spt is order, not expressed by time
%caution

e_a_list = [];
e_n_list = [];
del_e_list = [];

r_error_list= [];
energy_error_list = [];

an_dist_list = [];
nu_dist_list = [];

ext_total_list = [];
nu_total_list = [];

orbit_t_list = [];

interp_idx_list = [];
for spt= spt_list
    interp_idx = find(abs(interp_T-spt) < 0.0000001);
    interp_idx_list = [interp_idx_list interp_idx];
end

for interp_idx = interp_idx_list
    dist_error_list = [];
    
    m_ano = 2*pi/bi_period*interp_T(interp_idx);        %mean anomaly
    m_red = m_t{1}*m_t{2}/(m_t{1}+m_t{2});
    
%     %--------------------------center of mass------------------------------
%     x=0;y=0;z=0;mass=0;
%     vx=0;vy=0;vz=0;
%     for i = 1:length(pos_t)
%         x = x + m_t{i}*pos{i}(spt,1);
%         y = y + m_t{i}*pos{i}(spt,2);
%         z = z + m_t{i}*pos{i}(spt,3);
%         vx = vx + m_t{i}*vel{i}(spt,1);
%         vy = vy + m_t{i}*vel{i}(spt,2);
%         vz = vz + m_t{i}*vel{i}(spt,3);
%         mass = mass + m_t{i};
%     end
%     
%     x = x/mass;
%     y = y/mass;
%     z = z/mass;
%     vx = vx/mass;
%     vy = vz/mass;
%     vz = vz/mass;
%     
%     tmp_CM_pos = [x y z];
%     tmp_CM_vel = [vx vy vz];
%     
%     %in order to prove system is conserved 
%     CM_ang_m = sqrt(sum(tmp_CM_pos.^2))*m_red*sqrt(sum(tmp_CM_vel.^2));
%     
    %---------------------set up for eccentricity--------------------------
    %second object is moving to center of mass
    ang_m_a = dist(pos{1}(1,:),pos{2}(1,:))*m_t{2}*sqrt(sum(vel{2}(1,:).^2));
    ang_m_n = dist(interp_pos{1}(interp_idx,:),interp_pos{2}(interp_idx,:))*m_t{2}*sqrt(sum(interp_vel{2}(interp_idx,:).^2));
    TE_spt_a = E(1) + U(1);
    TE_spt_n = interp_E(interp_idx) + interp_U(interp_idx);
    
    alpha = G*m_t{1}*m_t{2};
    e_a = sqrt(1 + ((2*TE_spt_a*ang_m_a^2)/(m_red*alpha^2)));
    e_n = sqrt(1 + ((2*TE_spt_n*ang_m_n^2)/(m_red*alpha^2)));
    e_a_list = [e_a_list e_a];
    e_n_list = [e_n_list e_n];
    del_e_list = [del_e_list e_n-e_a];
  
    
    syms E_orbit
    equ1 = E_orbit -e_a*sin(E_orbit) -m_ano;
    E_orbit1 = vpasolve(equ1, E_orbit);

    syms theta 
    equ2 = cos(E_orbit1)*(1+e_a*cos(theta)) - e_a -cos(theta);
    theta1 = vpasolve(equ2, theta);
    
    an_dist = 1/(dist(pos_s,pos_p)*(1+e_a*cos(theta1)));
    nu_dist = dist(interp_pos{1}(interp_idx,:),interp_pos{2}(interp_idx,:));
    
    an_dist_list = [an_dist_list an_dist];
    nu_dist_list = [nu_dist_list nu_dist];
    
    a_error = an_dist - nu_dist;
    r_error = (an_dist - nu_dist)/an_dist; 

    r_error_list = [r_error_list r_error];
    
    %-----------------------------energy-----------------------------------
    %I gave initial total energy lower
    set_total = E(1) + U(1);
    ext_potential = - G*m_t{1}*m_t{2}/an_dist;
    ext_kinetic = set_total - ext_potential;
    
    ext_total = ext_kinetic + ext_potential;
    
    nu_total = interp_E(interp_idx)+interp_U(interp_idx);
    energy_error = (ext_total - nu_total)/ext_total;
    
    ext_total_list = [ext_total_list ext_total];
    nu_total_list = [nu_total_list nu_total];
    energy_error_list = [energy_error_list energy_error];

    orbit_t = interp_T(interp_idx);
    orbit_t_list = [orbit_t_list orbit_t];
end
 length(e_a_list)
 length(spt_list)

dT = diff(T);
eta(m) = max(dT)/min(dT)
eta_t  = dT/min(dT);

%eta as function of time
subplot(3,3,6);
    plot(T(1:end-1),eta_t);hold on;
    title('Eta as function of t');

%energy and error

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


figure('Name','primary energy','NumberTItle','off');
    subplot(2,3,1);
    plot(T,E1);hold on;
    title('primary kinetic');
    
    subplot(2,3,2);
    plot(T,U1);hold on;
    title('primary potential');
    
    subplot(2,3,3);
    plot(T,E1+U1);
    title('primary total energy');
    
    subplot(2,3,4);
    plot(T,E2);hold on;
    title('secondary star kinetic');
    
    subplot(2,3,5);
    plot(T,U2);hold on;
    title('secondary potential');
    
    subplot(2,3,6);
    plot(T,E2+U2);hold on;
    title('secondaty total energy');
 
% figure(4);
% plot(vel{1}(:,1),m_t{1}*vel{1}(:,1));
% xlabel('x');ylabel('y');
% title('P-Q diagram');
%-------------------------------3d plot------------------------------------
figure('Name','Motion of objects','NumberTitle','off');
    plot(pos{1}(:,1),pos{1}(:,2),'r');hold on;
    plot(pos{2}(:,1),pos{2}(:,2),'b');hold on;grid on;
    xlabel('x');ylabel('y');
    title('motion of component');

figure('Name','entire system"s energy','NumberTitle','off');
    subplot(3,1,1);
    plot(T,E);
    title('total kinetic');
    subplot(3,1,2);
    plot(T,U);
    title('total potential');
    subplot(3,1,3);
    plot(T,E+U);
    title('hamiltonian')

figure('Name','analytic and numerical eccentricity','NumberTitle','off');
    subplot(2,1,1);
    plot(orbit_t_list,e_a_list,'o');hold on;
    plot(orbit_t_list,e_n_list,'.','markersize',8);
    xlabel('orbits')
    ylabel('eccentricity')
    subplot(2,1,2);
    plot(orbit_t_list,del_e_list);
    xlabel('orbits')
    ylabel('discrepency of a and n eccentricity')
    
figure('Name','error','NumberTitle','off');
    subplot(2,1,1)
    plot(orbit_t_list,r_error_list)
    subplot(2,1,2)
    plot(orbit_t_list,energy_error_list)

figure('Name','exact and numerical solution','NumberTitle','off');
    subplot(2,1,1)
    plot(orbit_t_list,an_dist_list,'b');hold on;
    plot(orbit_t_list,nu_dist_list,'r');
    
    subplot(2,1,2)
    plot(orbit_t_list,ext_total_list,'b');hold on;
    plot(orbit_t_list,nu_total_list,'r');figure(123);
    
figure(111);
    subplot(3,1,1)
    plot(T,pos{1},'r');hold on;
    plot(interp_T,interp_pos{1},'o');
    subplot(3,1,2)
    plot(T,pos{1},'r');hold on;
    plot(interp_T,spline_pos{1},'o');
    subplot(3,1,3)
    plot(interp_T,(interp_pos{1}-spline_pos{1})/spline_pos{1});
    title('discrepancy of two methods')