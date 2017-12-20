function plot_2body_motion(pos,m_t,CM_mass)

figure('Name','Motion of objects','NumberTitle','off');
    plot(pos{1}(:,1),pos{1}(:,2),'r');hold on;
    plot(pos{2}(:,1),pos{2}(:,2),'b');grid on;
    
    rcm= (pos{1}*m_t{1} + pos{2}*m_t{2})/CM_mass;
    
    plot(rcm(:,1),rcm(:,2),'linewidth',4);
    plot(rcm(1,1),rcm(1,2),'o','markersize',14); 
    
    l=length(rcm(:,1));
    plot(rcm(l,1),rcm(l,2),'>','markersize',14);
    xlabel('x');ylabel('y');
    title('motion of component');
end
