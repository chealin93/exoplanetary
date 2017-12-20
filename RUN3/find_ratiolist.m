function [ratio_list,e_list,ea_list]=find_ratiolist()

G = 1;
m_s = [0.5];
m_p = [0.5];
m_t = {m_s m_p};

pos_S = [0 -0.5 0];
pos_P = -m_s/m_p*pos_S;
pos_t = {pos_S pos_P};

[CM_pos,CM_mass] = COM1(pos_t,m_t);
pos_s = pos_S - CM_pos;
pos_p = pos_P - CM_pos;

pos_t = {pos_s pos_p};

vel_sx = sqrt(G*m_t{2}^2/(m_t{1}+m_t{2})/dist(pos_P,pos_S));
vel_px = sqrt(G*m_t{1}^2/(m_t{1}+m_t{2})/dist(pos_S,pos_P));

e = 0;
ratio_list = [];
e_list = [];
ea_list = [];

while e < 0.91
    
    for ratio = 1:-0.0001:-2
        
        vel_S = [vel_sx 0 0];
        vel_P = ratio*[-vel_px 0 0];
        vel_t = {vel_S vel_P};
        
        CM_vel = COM2(vel_t,m_t);

        vel_s = vel_S - CM_vel;
        vel_p = vel_P - CM_vel;
        vel_t = {vel_s vel_p};

        fir_ang = dist(CM_pos,pos_t{2}(1,:))*sqrt(sum(vel_t{2}(1,:).^2));
        fir_E = 1/2*m_t{1}*sum(vel_t{1}(1,:).^2)...
            + 1/2*m_t{2}*sum(vel_t{2}(1,:).^2);
        fir_U = - G*m_t{1}*m_t{2}/dist(pos_s,pos_p);
        fir_H = fir_E + fir_U;
        alpha = m_t{1}*m_t{2};
        two_m_red = m_t{1}*m_t{2}/(m_t{1}+m_t{2});
        e_a = sqrt(1 + ((2*fir_H*fir_ang^2)/(two_m_red*alpha^2)));
        
       if abs(abs(e_a) - e) < 0.0001
           ratio_list = [ratio_list ratio];
           ea_list = [ea_list e_a];
           break
       end
    end
    e_list = [e_list e];
    e = e + 0.02;
end
end



            
        