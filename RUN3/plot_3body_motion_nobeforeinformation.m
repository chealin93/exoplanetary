function plot_3body_motion(pos,m_t)

figure('Name','Motion of objects','NumberTitle','off');
    subplot(2,2,1)
    plot(pos{1}(1:end/100,1),pos{1}(1:end/100,2),'r');hold on;
    plot(pos{2}(1:end/100,1),pos{2}(1:end/100,2),'b');grid on;
    plot(pos{3}(1:end/100,1),pos{3}(1:end/100,2),'k');
    
    subplot(2,2,2)
    plot(pos{1}(end/100:end*2/100,1),pos{1}(end/100:end*2/100,2),'r');hold on;
    plot(pos{2}(end/100:end*2/100,1),pos{2}(end/100:end*2/100,2),'b');grid on;
    plot(pos{3}(end/100:end*2/100,1),pos{3}(end/100:end*2/100,2),'k');
    
    subplot(2,2,3)
    plot(pos{1}(end*2/100:end*3/100,1),pos{1}(end*2/100:end*3/100,2),'r');hold on;
    plot(pos{2}(end*2/100:end*3/100,1),pos{2}(end*2/100:end*3/100,2),'b');grid on;
    plot(pos{3}(end*2/100:end*3/100,1),pos{3}(end*2/100:end*3/100,2),'k');
    
    subplot(2,2,4)
    plot(pos{1}(end*3/100:end*4/100,1),pos{1}(end*3/100:end*4/100,2),'r');hold on;
    plot(pos{2}(end*3/100:end*4/100,1),pos{2}(end*3/100:end*4/100,2),'b');grid on;
    plot(pos{3}(end*3/100:end*4/100,1),pos{3}(end*3/100:end*4/100,2),'k');
    
    rcm = (pos{1}*m_t{1} + pos{2}*m_t{2} +pos{3}*m_t{3})/(m_t{1}+m_t{2}+m_t{3});
    
    plot(rcm(1:100,1),rcm(1:100,2),'linewidth',4);
    plot(rcm(1,1),rcm(1,2),'o','markersize',14); 
    
    l=length(rcm(1:100,1));
    plot(rcm(l,1),rcm(l,2),'>','markersize',14);
    xlabel('x');ylabel('y');
    title('motion of component');
end
