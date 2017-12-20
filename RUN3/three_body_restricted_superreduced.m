%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@                          exoplanetary system                           @ 
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

close all;
format long;
clear all;

global eps;
global m_t;
global number_of_planet;
global tot_num;
global G;
global R;
global omega;
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
separation2 = 0.8;   
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
m_ratio = 1;
m_p = m_s*m_ratio;
m_t = {m_s m_p};         %normalized by earth mass

pos_S = [0 -0.5 0];
pos_P = -m_s/m_p*pos_S;
pos_t = {pos_S pos_P};

[CM_pos,CM_mass] = COM1(pos_t,m_t);

pos_s = pos_S - CM_pos;
pos_p = pos_P - CM_pos;

vel_sx = sqrt(G*m_t{2}^2/(m_t{1}+m_t{2})/dist(pos_P,pos_S));
vel_px = sqrt(G*m_t{1}^2/(m_t{1}+m_t{2})/dist(pos_S,pos_P));
        
vel_S = 100*[vel_sx 0 0];
vel_P = 100*[-vel_px 0 0];
vel_t = {vel_S vel_P};

CM_vel = COM2(vel_t,m_t);

vel_s = vel_S - CM_vel;
vel_p = vel_P - CM_vel;

%-------------------------------------------------------------------
m_c = [1e-8];

phi = 0;
dx2 = separation2*sind(phi);
dy2 = separation2*cosd(phi);
        
pos_vec = [dx2 dy2 0];
        
x2 = CM_pos(1)+ dx2;
y2 = CM_pos(2)+ dy2;
pos_C = [x2 y2 0];

mag_pos = sqrt(sum(pos_vec.^2));
vel_perpen_unit = [dy2 -dx2 0]/mag_pos;
pos_c = pos_C - CM_pos;


vel_cx = sqrt(G*CM_mass/dist(CM_pos,pos_C));
vel_C = -vel_cx*vel_perpen_unit;
vel_c = vel_C - CM_vel;
            
pos_t = {pos_s pos_p pos_c};
vel_t = {vel_s vel_p vel_c};
m_t = {m_s m_p m_c};
CM_pos = CM_pos - CM_pos;

%------for eccentricity---------

fir_ang = dist(CM_pos,pos_t{2}(1,:))*sqrt(sum(vel_t{2}(1,:).^2));
%fir_ang = specific angular momentum
fir_E = 1/2*m_t{1}*sum(vel_t{1}(1,:).^2)...
    + 1/2*m_t{2}*sum(vel_t{2}(1,:).^2);
fir_U = - G*m_t{1}*m_t{2}/dist(pos_s,pos_p);
fir_H = fir_E + fir_U;
alpha = m_t{1}*m_t{2};
two_m_red = m_t{1}*m_t{2}/(m_t{1}+m_t{2});
e_a = sqrt(1 + ((2*fir_H*fir_ang^2)/(two_m_red*alpha^2)))

a = dist(CM_pos,pos_p);
j = abs(sqrt(G*CM_mass*a*(1-e_a^2)));
A = abs(pi*a^2*sqrt(1-e_a^2));
pperiod = 2*A/j

%---------global parameter in ode45 solver--------
omega = 2*pi/pperiod;
R = dist(CM_pos,pos_p);

%--------------------------------------------------------------------------
% integration by using ode45 solver
%--------------------------------------------------------------------------
m=1; eps=0; t0 =0; tf=pperiod*80;
tot_num = length(pos_t); n = tot_num;

y0 = [pos_t{3} vel_t{3}];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T,Y] = ode45('several_planet_moon_rk_superreduced',[t0,tf],y0,options);

%--------------------------------------------------------------------------
% separating object's position and velocity to cell easily
% additionaly,distance of all objects 
%--------------------------------------------------------------------------

pos = {};
vel = {};
pos{3} = Y(:,1:3);
vel{3} = Y(:,1:3);

