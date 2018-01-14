%--------------------------------------
% Center of mass 
% 1. positon 2. total mass of system
%--------------------------------------

function  [CM_pos,CM_mass] = COM1(pos_t,m_t)

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
CM_mass = tot_mass;
end

