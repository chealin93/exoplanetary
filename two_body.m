%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@                          exoplanetary system                           @ 
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

close all;
clear all;
format shortEng;

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
separation2 = 1.5;   
inclination1 = 0;       
inclination2 = 0;       

%-----------------------------set exosystem--------------------------------
%recognize the exoplanetary system and properties based on solar system
%proportional expression (standard 'the erath')

%--------------------------------------------------------------------------
% planetary system on the center of mass frame
%--------------------------------------------------------------------------
G = 1;
n_p = number_of_planet + 2;
m_s = [0.5];
m_p = [0.5];
m_t = {m_s m_p};         %normalized by earth mass

se_t = {separation1};
inclination_t = {inclination1};

pos_S = [0 0 0];
dy1 = cosd(inclination1)*se_t{1};
dz1 = sind(inclination1)*se_t{1};
pos_P = [0 dy1 dz1];
pos_t = {pos_S pos_P};

[CM_pos,CM_mass] = COM1(pos_t,m_t);

pos_s = pos_S - CM_pos;
pos_p = pos_P - CM_pos;

rs = sqrt(sum(pos_s.^2));
rp = sqrt(sum(pos_p.^2));

vel_sx = sqrt(G*m_t{2}*rs/(dist(pos_P,pos_S)));
vel_px = sqrt(G*m_t{1}*rp/(dist(pos_S,pos_P)));

vel_S = [vel_sx 0 0];
vel_P = [-vel_px 0 0];
vel_t = {vel_S vel_P};

CM_vel = COM2(vel_t,m_t);

vel_s = vel_S - CM_vel;
vel_p = vel_P - CM_vel;

pos_t = {pos_s pos_p};
vel_t = {vel_s vel_p};
CM_pos = CM_pos - CM_pos;

fir_ang = dist(CM_pos,pos_t{2}(1,:))*sqrt(sum(vel_t{2}(1,:).^2));
fir_E = 1/2*m_t{1}*sum(vel_t{1}(1,:).^2)...
    + 1/2*m_t{2}*sum(vel_t{2}(1,:).^2);
fir_U = - G*m_t{1}*m_t{2}/dist(pos_s,pos_p);
fir_H = fir_E + fir_U;
alpha = m_t{1}*m_t{2};
two_m_red = m_t{1}*m_t{2}/(m_t{1}+m_t{2});
e_a = sqrt(1 + ((2*fir_H*fir_ang^2)/(two_m_red*alpha^2)))


%--------------------------------------------------------------------------
% integration by using ode45 solver
%--------------------------------------------------------------------------
m=1; eps=0; t0 =0; tf=100;
tot_num = length(pos_t); n = tot_num;

y0 = [pos_t{1:n} vel_t{1:n}];length(y0);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T,Y] = ode45('several_planet_moon_rk',[t0,tf],y0,options);

%--------------------------------------------------------------------------
% separating object's position and velocity to cell easily
% additionaly,distance of all objects 
%--------------------------------------------------------------------------

k=1;j=3*n+1;
pos={};vel={};
for i = 1:n
    pos{i}=Y(:,k:k+2);
    k=k+3;
    vel{i}=Y(:,j:j+2);
    j=j+3;
end

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

%--------------------------------------------------------------------------
% energy calculation
%--------------------------------------------------------------------------
[all_E,all_U,tot_E,tot_U] = energyc(vel,m_t,D_t,R_m_t,T);
H = tot_E + tot_U;

%--------------------------------------------------------------------------
% determining orbiting period 
%--------------------------------------------------------------------------
% from center of mass to position2 object
[semi_a,period] = kepler3rd(m_t,H)

spt = period;
spt_list = [];
while tf > spt
    spt_list = [spt_list spt];
    spt = spt + period;
end
period_list = spt_list;
period_point_number = length(period_list);

spt_idx_list = [];
for spt = spt_list
    [dif idx] = min(abs(T-spt));
    spt_idx_list = [spt_idx_list idx];
end


%-------------plot energy and numerical result-----------------
plot_analysis(T,pos,vel,m_t,D_t,tot_E,tot_U,m)   
plot_2body_motion(pos,m_t,CM_mass)

%--------------------------------------------------------------------------
% exact solution
%--------------------------------------------------------------------------
%spt is specific time, when you want to give the specific time,
%you must call like T(spt), spt is order, not expressed by time
%caution
e_a_list = [];
e_n_list = [];
del_e_list = [];

ang_a_list = [];
ang_n_list = [];
ang_e_list = [];
energy_error_list = [];

total_a_list = [];
total_n_list = [];

orbit_list = [];
spt_list = [];
orbit = 1;

spt_idx = [];
[~,c]=ismember(spt_list,T);

H = tot_E + tot_U;

for spt = spt_idx_list
    
    m_red = m_t{1}*m_t{2}/CM_mass; 
    %---------------------set up for eccentricity--------------------------
    %second object is moving to center of mass
    ang_m_a = dist(CM_pos,pos{2}(1,:))*sqrt(sum(vel{2}(1,:).^2));
    ang_m_n = dist(CM_pos,pos{2}(spt,:))*sqrt(sum(vel{2}(spt,:).^2));
    
    ang_a_list = [ang_a_list ang_m_a];
    ang_n_list = [ang_n_list ang_m_n];
    ang_error = ang_m_a - ang_m_n;
    ang_e_list = [ang_e_list ang_error];
    
    TE_spt_a = H(1);
    TE_spt_n = tot_E(spt) + tot_U(spt);
    
    alpha = G*m_t{1}*m_t{2};
    e_a = sqrt(1 + ((2*TE_spt_a*ang_m_a^2)/(m_red*alpha^2)))
    e_n = sqrt(1 + ((2*TE_spt_n*ang_m_n^2)/(m_red*alpha^2)));
    
    e_a_list = [e_a_list e_a];
    e_n_list = [e_n_list e_n];
    del_e =  ang_m_a - ang_m_n;
    del_e_list = [del_e_list del_e];
  
    %-----------------------------energy-----------------------------------
    %I gave initial total energy lower
    ext_kinetic = 1/2*m_t{1}*sum(vel{1}(1,:).^2) + 1/2*m_t{2}*sum(vel{2}(1,:).^2);
    ext_potential = - G*m_t{1}*m_t{2}/dist(pos{1}(1,:),pos{2}(1,:));
    ext_total = ext_kinetic + ext_potential;
    
    nu_total = tot_E(spt)+tot_U(spt);
    energy_error = ext_total - nu_total;
    
    total_a_list = [total_a_list ext_total];
    total_n_list = [total_n_list nu_total];
    energy_error_list = [energy_error_list energy_error];

    orbit = orbit+1;
    spt_list = [spt_list spt];
    orbit_list = [orbit_list orbit];
end


 %------------------------------------------------------------------------
 % plotting region
 %-------------------------------------------------------------------------

plot_energy(T,all_E,all_U)  
plot_H(T,tot_E,tot_U)
plot_e_E1(orbit_list,e_a_list,ang_a_list,ang_n_list,ang_e_list,e_n_list,del_e_list,total_a_list,total_n_list,energy_error_list)
plot_e_E(orbit_list,e_a_list,ang_a_list,ang_n_list,ang_e_list,e_n_list,del_e_list,total_a_list,total_n_list,energy_error_list)