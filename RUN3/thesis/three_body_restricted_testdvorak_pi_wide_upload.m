%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@                          exoplanetary system                           @ 
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
% This code is exoplanetary N-body simulator for detection possibility of
% various type of exoplanetary systems.
%
% Simulator is made using MATLABR2016b, you want to make your style.

format long;
clear all;
close all;

global eps;
global m_t;
global number_of_planet;
global tot_num;

%--------------------set exoplanetary system-------------------------------
% elliptic restricted three-body problem code
% planetary system on the center of mass frame
% distance, mass is normalized by AU, solar mass
% x,y,z coordinate system
% always equal massive binary
% mass of the third body is up to 10^-4 from 10^-8

G = 1;
n_p = number_of_planet + 2;
m_s = [0.5];
m_ratio = 1;
m_p = m_s*m_ratio;
m_t = {m_s m_p};        

pos_S = [0 -0.5 0];
pos_P = -m_s/m_p*pos_S;
pos_t = {pos_S pos_P};

[CM_pos,CM_mass] = COM1(pos_t,m_t); 
tmp_CM = COM1(pos_t,m_t);

pos_s = pos_S - CM_pos;
pos_p = pos_P - CM_pos;

vel_sx = sqrt(G*m_t{2}^2/(m_t{1}+m_t{2})/dist(pos_P,pos_S));
vel_px = sqrt(G*m_t{1}^2/(m_t{1}+m_t{2})/dist(pos_S,pos_P));

%==============================================
% eccentricity is changed by varying velocity
%==============================================

[ratio_list,b,c] = find_ratiolist();
%ratio_list is list for finding eccentricity from 0 to 0.9 
%with step size 0.02
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

order = 1;
drangei = 0.8; % initial minimum distance
drangef = 3.5; % initial maximum distance

for ratio = ratio_list
    vel_S = [vel_sx 0 0]
    vel_P = ratio*[-vel_px 0 0]
    vel_t = {vel_S vel_P};

    CM_vel = COM2(vel_t,m_t);

    vel_s = vel_S - CM_vel;
    vel_p = vel_P - CM_vel;
    %----------------------------------------------------------------------
    m_c = [1e-8]; % third body mass
    
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
    
    if rem(order,6) == 0
        drangef = drangef + 1;
    end
    
    if order >= 24 && rem(order,6) == 0
        drangei = drangei + 1;
    end
    
    %separation2 is initial distance between the third body and center of
    %mass of two primaries
    
    for separation2 = drangei:0.05:drangef
        in = 0;
        out = 0;
        in_ang_list = [];
        in_dist_list = [];
        out_ang_list = [];
        out_dist_list = [];
        
        % phi is angle to the line connecting with two primaries
        
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
            % initial direction of velocity of the third body is
            % perpendicular to the center of mass of two primaries

            pos_t = {pos_s pos_p pos_c};
            vel_t = {vel_s vel_p vel_c};
            m_t = {m_s m_p m_c};
        
            CM_pos = CM_pos - CM_pos;
        
            %------------------finding eccentricitiy----------------------- 
            fir_ang = dist(CM_pos,pos_t{2}(1,:))*sqrt(sum(vel_t{2}(1,:).^2));
            fir_E = 1/2*m_t{1}*sum(vel_t{1}(1,:).^2)...
                + 1/2*m_t{2}*sum(vel_t{2}(1,:).^2);
            fir_U = - G*m_t{1}*m_t{2}/dist(pos_s,pos_p);
            fir_H = fir_E + fir_U;
            alpha = m_t{1}*m_t{2};
            two_m_red = m_t{1}*m_t{2}/(m_t{1}+m_t{2});
            e_a = sqrt(1 + ((2*fir_H*fir_ang^2)/(two_m_red*alpha^2)))
            % e_a is numerical eccentricity from energy 
            
            pperiod = 2*pi*dist(CM_pos,pos_p)/sqrt(sum(vel_P.^2));
            omega = 2*pi/pperiod;
            R = dist(CM_pos,pos_p);
            %--------------------------------------------------------------
            % integration by using ode45 solver
            %--------------------------------------------------------------
            % at least, integration time is 500 times period of one of the
            % two primaries
            
            m=1; eps=0; t0 =0; tf=pperiod*500;
            tot_num = length(pos_t); n = tot_num;

            y0 = [pos_t{1:n} vel_t{1:n}];length(y0);
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            % you can choose numerical error up to 10^-12 that is adaptive
            % ode 45 solver 
            [T,Y] = ode45('several_planet_moon_rk',[t0,tf],y0,options);
            
            % T is time list
            % Y is integrated position and velocity

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
        
            % after 500 times period integration, we clasify unstable
            % orbita and stable orbit by judging distance between the third
            % body and center or mass
            
            if dist(pos{3}(end,:),CM_pos) > 30
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
    
    % after previous process, choose maximum value in UCO_init_dist_list
    % and minimum value in LCO_init_dist_list
    
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
    order = order + 1;
end

