%---------------------------------------------------------
% given eccentricity, calculating third initial velocity
%---------------------------------------------------------

function [given_vel1,given_vel2] = given_vel(given_e,m_t,pos_t,vel_t)

G = 1;
m_red = m_t{1}*m_t{2}*m_t{3}/(m_t{1}+m_t{2}+m_t{3});
alpha = G*m_t{1}*m_t{2}*m_t{3};


x=0;y=0;z=0;
tot_mass=0;

for i=1:length(pos_t)
    tot_mass = tot_mass + m_t{i};
    x = x + m_t{i}*pos_t{i}(1);
    y = y + m_t{i}*pos_t{i}(2);
    z = z + m_t{i}*pos_t{i}(3);
end

x = x/tot_mass;
y = y/tot_mass;
z = z/tot_mass;

CM_pos = [x y z];

dist_com_2 = dist(CM_pos,pos_t{2}(1,:));

ext_E = 1/2*m_t{1}*sum(vel_t{1}(1,:).^2) + 1/2*m_t{3}*sum(vel_t{3}(1,:).^2);
ext_U = - G*m_t{1}*m_t{2}/dist(pos_t{1},pos_t{2}) - G*m_t{2}*m_t{3}/dist(pos_t{2},pos_t{3}) - G*m_t{3}*m_t{1}/dist(pos_t{3},pos_t{1}); 
tot_E = ext_E + ext_U;

syms T
equ1 = m_t{2}^3*dist_com_2^2*T^2 + 2*m_t{2}^2*dist_com_2^2*tot_E*T - m_red*alpha*(given_e^2 -1);
vv = solve(equ1,T);
given_vel = sqrt(vv);
given_vel1 = double(subs(given_vel(1)));
given_vel2 = double(subs(given_vel(2)));
end
