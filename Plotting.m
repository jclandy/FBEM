function Plotting(topo_plot,echo_plot,PosT,t,P_t_ml,P_t_full,P_t_ml_comp,N_b,epsilon_b)
    
    if topo_plot>0
        
        figure(1); clf(1); 
        TRI = delaunay(PosT(:,1),PosT(:,2));
        trisurf(TRI,PosT(:,1),PosT(:,2),PosT(:,3),'edgecolor','none')
        colorbar('southoutside')
        view(-45,80)
        xlim([-1000 1000])
        ylim([-1000 1000])
        title('Snow-covered sea ice surface tetrahedral mesh')
        xlabel('X [m]')
        ylabel('Y [m]')
        zlabel('Z [m]')
        drawnow;
        
    else
    end
    
    if echo_plot>0
        
        figure(2); clf(2);
        subplot(2,2,1)
        plot(t*1e9,P_t_ml/(max(P_t_ml)),'linewidth',2)
        xlim([-20 80])
        ylim([0 1])
        grid on
        title('Multi-looked power echo')
        xlabel('Time [ns]')
        ylabel('Normalized Power')
        
        subplot(2,2,2)
        plot(t*1e9,P_t_full(:,1:4:N_b),'linewidth',1)
        xlim([-20 80])
        grid on
        title('Single-look power echoes')
        xlabel('Time [ns]')
        ylabel('Power [W]')
        legend(sprintfc('%d',1:4:N_b),'Location','eastoutside')
        
        subplot(2,2,3)
        m = -(N_b-1)/2:(N_b-1)/2;
        imagesc(t*1e9,m*epsilon_b*180/pi,P_t_full')
        xlim([-20 80])
        title('Delay-Doppler Map')
        xlabel('Time [ns]')
        ylabel(['Look Angle [' char(176) ']'])
        colorbar
        
        subplot(2,2,4)
        plot(t*1e9,P_t_ml_comp,'linewidth',1)
        xlim([-20 80])
        grid on
        title('Multi-looked component echoes')
        xlabel('Time [ns]')
        ylabel('Power [W]')
        legend('Snow Surf','Snow Vol','Ice Surf','Lead/Pond Surf','Location','northeast')
        drawnow;
        
    else
    end

end
