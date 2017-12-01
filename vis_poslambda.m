function = vis_poslambda(num_backsteps,spikes,varargin)

    figure(); clf; hold on;
    set(gcf,'units','points','position',[100,100,1000,400])

    subplot(1,2,1); hold on;
    y_new = flipud(y_new);
    x_new = fliplr(x_new);

    lambda3(find(x_new.^2 + y_new.^2 > 1))=nan;  % what is this line doing?

    % plot lambda as a function position over this grid
    h_mesh = mesh(x_new,y_new,lambda3,'AlphaData',0);
    get(h_mesh,'AlphaData');
    set(h_mesh,'AlphaData',0);
    hold on;
    plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
    xlabel('x position [m]'); ylabel('y position [m]');