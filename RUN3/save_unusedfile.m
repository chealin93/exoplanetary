
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