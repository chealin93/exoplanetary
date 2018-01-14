
function dy = f(t,y)

global eps;
global m_t;
global tot_num;
global omega;
global G;

n = tot_num;
p = 3;
phi = omega*t;

X=y(1:9);
V=y(10:18);

vel_t = [];
acc_t = [];
dy = [];
%--------position of n-body----------

pos_s = X(1:3)';
pos_p = X(4:6)';
pos_c = X(7:9)';
pos_t = {pos_s pos_p pos_c};
%----------velocity of n-body---------

[CM_pos,CM_mass]=COM1(pos_t(1:2),m_t(1:2));

dx1 = pos_s(1);
dy1 = pos_s(2);

pos_vec = [dx1 dy1 0];           
mag_pos = sqrt(sum(pos_vec.^2));
vel_perpen_unit = [dy1 -dx1 0]/mag_pos;
vel_sx = sqrt(G*m_t{1}/dist(pos_s,pos_p)*dist(CM_pos,pos_s));
vel_s = -vel_sx*vel_perpen_unit;
            
dx2 = pos_p(1);
dy2 = pos_p(2);
pos_vec = [dx2 dy2 0];
vel_px = sqrt(G*m_t{2}/dist(pos_p,pos_s)*dist(CM_pos,pos_p));
mag_pos = sqrt(sum(pos_vec.^2));
vel_perpen_unit = [dy2 -dx2 0]/mag_pos;
vel_p = -vel_px*vel_perpen_unit;

vel_c = V(7:9)';
vel_t = [vel_s vel_p vel_c];
%----acceleration of thirdbody-----
a1 = - m_t{2}*G*(pos_t{1}-pos_t{2})/(sqrt(sum((pos_t{1}-pos_t{2}).^2))+eps).^p;
a2= - a1;
a3 = - m_t{1}*(pos_t{3}-pos_t{1})/(sqrt(sum((pos_t{3}-pos_t{1}).^2))+eps).^p...
    - m_t{2}*(pos_t{3}-pos_t{2})/(sqrt(sum((pos_t{3}-pos_t{2}).^2))+eps).^p;
acc_t = [a1 a2 a3];

%----------making dy---------------

dy = [vel_t acc_t];
dy = dy';
end