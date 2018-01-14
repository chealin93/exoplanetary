%---------------------------------------------
% center of mass velocity by using momentum
%---------------------------------------------

function CM_vel = COM2(vel_t,m_t)

tot_mass = 0;
vx = 0; vy = 0; vz = 0;
for i = 1:length(vel_t)
    vx = vx + m_t{i}*vel_t{i}(1);
    vy = vy + m_t{i}*vel_t{i}(2);
    vz = vz + m_t{i}*vel_t{i}(3);
    tot_mass = tot_mass + m_t{i};
end

vx = vx/tot_mass;
vy = vy/tot_mass;
vz = vz/tot_mass;

CM_vel = [vx vy vz];
end


