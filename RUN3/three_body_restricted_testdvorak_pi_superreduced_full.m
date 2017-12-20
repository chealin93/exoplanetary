%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@                          exoplanetary system                           @ 
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
format long;
clear all;
close all;

global eps;
global m_t;
global number_of_planet;
global tot_num;
global omega;
global R;
global G;
global separation2
%-------------------------input parameter----------------------------------
%input function
%when you put a value in each input function, you must enter exact value
%not type space!!!!!

number_of_planet = 0; %input('put the number of additional planets');
yesno = 2; %input('do you want to generate random number of moons? answer yes(1) or no(2)');

if yesno == 2
   input_num_m = 0; %input('put the additional moon to each planets');
end

separation1 = 1; 
separation2 = 0;
inclination1 = 0;       
inclination2 = 0;       

%-----------------------------set exosystem--------------------------------
%recognize the exoplanetary system and properties based on solar system
%proportional expression (standard 'the erath')

%================================================
% planetary system on the center of mass frame
%================================================
G = 1;
n_p = number_of_planet + 2;
m_s = [0.5];
m_ratio = 1;
m_p = m_s*m_ratio;
m_t = {m_s m_p};         %normalized by earth mass

se_t = {separation1 separation2};
inclination_t = {inclination1 inclination2};

pos_S = [0 -0.5 0];
pos_P = -m_s/m_p*pos_S;
pos_t = {pos_S pos_P};

[CM_pos,CM_mass] = COM1(pos_t,m_t);
tmp_CM = COM1(pos_t,m_t);

pos_s = pos_S - CM_pos;
pos_p = pos_P - CM_pos;

vel_sx = sqrt(G*m_t{1}/dist(pos_S,pos_P)*dist(CM_pos,pos_S));
vel_px = sqrt(G*m_t{2}/dist(pos_P,pos_S)*dist(CM_pos,pos_P));

%==============================================
% eccentricity is changed by varying velocity
%==============================================

[ratio_list] = find_ratiolist();
angle_list = [0 45 90 135 180];

tot_in_ang_cell={};
tot_in_dist_cell={};
tot_out_ang_cell={};
tot_out_dist_cell={};

UCO_val = [];
LCO_val = [];
eccen_UCO = [];
eccen_LCO = [];

h = 1;
hh = 1;

limit_eccentricity = [];

vel_list = [];
vel_UCO = [];
vel_LCO = [];
for velratio = 1:10:30
    vel_S = velratio*[vel_sx 0 0];
    vel_P = velratio*[-vel_px 0 0];
    vel_t = {vel_S vel_P};

    CM_vel = COM2(vel_t,m_t);

    vel_s = vel_S - CM_vel;
    vel_p = vel_P - CM_vel;
    %--------------------------------------------------------------------------
    m_c = [1e-8];
    
    in_ang_cell = {};
    in_dist_cell = {};
    out_ang_cell = {};
    out_dist_cell = {};

    LCO_init_pos_list = [];
    LCO_init_ang_list = [];
    LCO_init_dist_list = [];
                
    UCO_init_pos_list = [];
    UCO_init_ang_list = [];
    UCO_init_dist_list = [];

    p = 1;
    for separation2 = 0.8:0.05:3
        in = 0;
        out = 0;
        in_ang_list = [];
        in_dist_list = [];
        out_ang_list = [];
        out_dist_list = [];
    
        for phi = angle_list
        
            dx2 = separation2*sind(phi);
            dy2 = separation2*cosd(phi);
        
            pos_vec = [dx2 dy2 0];
        
            x2 = CM_pos(1)+ dx2;
            y2 = CM_pos(2)+ dy2;
        
            pos_C = [x2 y2 0];
        
            mag_pos = sqrt(sum(pos_vec.^2));
            vel_perpen_unit = [dy2 -dx2 0]/mag_pos;
            pos_c = pos_C - CM_pos;

            rc = sqrt(sum(pos_c.^2));
            vel_cx = sqrt(G*CM_mass/dist(CM_pos,pos_C));
            vel_C = -vel_cx*vel_perpen_unit;
            vel_c = vel_C - CM_vel;

            pos_t = {pos_s pos_p pos_c};
            vel_t = {vel_s vel_p vel_c};
            m_t = {m_s m_p m_c};
        
            CM_pos = CM_pos - CM_pos;
            
            fir_ang = dist(CM_pos,pos_t{2}(1,:))*sqrt(sum(vel_t{2}(1,:).^2));
            fir_E = 1/2*m_t{1}*sum(vel_t{1}(1,:).^2)...
                + 1/2*m_t{2}*sum(vel_t{2}(1,:).^2);
            fir_U = - G*m_t{1}*m_t{2}/dist(pos_s,pos_p);
            fir_H = fir_E + fir_U;
            alpha = m_t{1}*m_t{2};
            two_m_red = m_t{1}*m_t{2}/(m_t{1}+m_t{2});
            e_a = sqrt(1 + ((2*fir_H*fir_ang^2)/(two_m_red*alpha^2)))
            
            pperiod = 2*pi*dist(CM_pos,pos_p)/sqrt(sum(vel_P.^2));

            %---------global parameter in ode45 solver--------
            omega = 2*pi/pperiod
            R = dist(CM_pos,pos_p);
        %--------------------------------------------------------------------------
        % integration by using ode45 solver
        %--------------------------------------------------------------------------
           m=1; eps=0; t0 =0; tf=pperiod*500;
            tot_num = length(pos_t); n = tot_num;

            G = m_t{2}*sum(vel_P.^2)/dist(CM_pos,pos_p)*dist(pos_s,pos_p)
            
            y0 = [pos_t{1:n} vel_t{1:n}];length(y0);
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            [T,Y] = ode45('Copy_of_several_planet_moon_rk_superreduced',[t0,tf],y0,options);

        %--------------------------------------------------------------------------
        % separating object's position and velocity to cell easily
        % additionaly,distance of all objects 
        %--------------------------------------------------------------------------
            
            pos={};vel={};
            pos{1}=Y(:,1:3);
            vel{1}=Y(:,10:12);
            
            pos{2}=Y(:,4:6);
            vel{2}=Y(:,13:15);
            
            pos{3}=Y(:,7:9);
            vel{3}=Y(:,16:18);
            
            l=length(T)
    D=[];D_t={};q=1;
    for i = 1:2
        for j = i+1:3
            for k =1:l
                D(k) = sqrt(sum((pos{i}(k,:)-pos{j}(k,:)).^2));
            end
                
            D_t{q} = D;
            R_m_t{q} = m_t{i}*m_t{j}; 
            q = q+1;
        end
    end

            [all_E,all_U,tot_E,tot_U] = energyc(vel,m_t,D_t,R_m_t,T);
            plot_energy(T,all_E,all_U)
            
            if dist(pos{1}(end,:),CM_pos) > 30
                out_ang_list = [out_ang_list phi];
                out_dist_list = [out_dist_list separation2];
                out = out + 1;
            else
                in_ang_list = [in_ang_list phi];
                in_dist_list = [in_dist_list separation2];
                in = in + 1;
            end 
        end
        
        in_ang_cell{p} = in_ang_list;
        in_dist_cell{p} = in_dist_list;
        out_ang_cell{p} = out_ang_list;
        out_dist_cell{p} = out_dist_list;
        p = p+1;
    
        if out == 5
            LCO_init_pos_list = [LCO_init_pos_list pos_c];
            LCO_init_ang_list = [LCO_init_ang_list phi];
            LCO_init_dist_list = [LCO_init_dist_list separation2];
        end
        
        if in == 5 
            UCO_init_pos_list = [UCO_init_pos_list pos_c];
            UCO_init_ang_list = [UCO_init_ang_list phi];
            UCO_init_dist_list = [UCO_init_dist_list separation2];
        end
    end
    
    if length(UCO_init_dist_list) >= 1
        UCO_val = [UCO_val min(UCO_init_dist_list)];
        tot_in_ang_cell{h} = in_ang_cell;
        tot_in_dist_cell{h} = in_dist_cell;
        eccen_UCO = [eccen_UCO e_a];
        h = h+1;
    end
    
    if length(LCO_init_dist_list) >= 1
        LCO_val = [LCO_val max(LCO_init_dist_list)];
        tot_out_ang_cell{h} = out_ang_cell;
        tot_out_dist_cell{h} = out_dist_cell;
        eccen_LCO = [eccen_LCO e_a];
        hh = hh + 1;
    end
    vel_UCO = [vel_UCO UCO_val];
    vel_LCO = [vel_LCO LCO_val];
    vel_list = [vel_list velratio];
end

figure(2);
plot(pos{1}(:,1),pos{1}(:,2));hold on;
plot(pos{2}(:,1),pos{2}(:,2));
plot(pos{3}(:,1),pos{3}(:,2));
xlabel('x')
ylabel('y')