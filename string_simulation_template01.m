function string_simulation_template01()

    close all
    clc

    % String parameters
    num_masses = 3;
    total_mass = 1; % kg
    tension_force = 10; % N
    string_length = 1; % m
    damping_coeff = 0.01; % Ns

    dx = string_length/(num_masses + 1);

    amplitude_Uf = 0.005; % m
    omega_Uf = 8.39; % rad/s

    Uf_func = @(t) amplitude_Uf*cos(omega_Uf*t);
    dUfdt_func = @(t) -amplitude_Uf*omega_Uf*sin(omega_Uf*t);

    % Parameters for the rate function
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;

    my_rate_func = @(t,V) string_rate_func01(t,V,string_params);

    % Initial pos and velo
    xlist = linspace(0,string_length,num_masses + 2);
    U0 = amplitude_Uf*xlist(2:end-1)'/string_length;
    dUdt0 = zeros(num_masses,1);
    V0 = [U0; dUdt0];

    dt = 1/600;
    tlist = 0:dt:12;

    Vlist = zeros(2*num_masses,length(tlist));
    Vlist(:,1) = V0;

    % RK Integration
    for i = 1:length(tlist)-1
        t = tlist(i);
        V = Vlist(:,i);

        k1 = my_rate_func(t,V);
        k2 = my_rate_func(t + dt/2,V + dt*k1/2);
        k3 = my_rate_func(t + dt/2,V + dt*k2/2);
        k4 = my_rate_func(t + dt,V + dt*k3);

        Vlist(:,i+1) = V + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % fixed and moving endpoints
    Ulist = [zeros(1,length(tlist));
             Vlist(1:num_masses,:);
             Uf_func(tlist)];

    % animation setup
    max_displacement = max(abs(Ulist(:)));
    y_limit = max(0.05,1.2*max_displacement);

    fig = figure('Color','w','Units','pixels', 'Position',[100 100 1440 1080]);

    string_plot = plot(xlist,Ulist(:,1),'-o','LineWidth',2,'MarkerSize',9,'Color',[0 0.35 0.7],'MarkerFaceColor',[0 0.35 0.7],'MarkerEdgeColor',[0 0.35 0.7]);

    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 20;

    title('Vibration of a String with Discrete Masses','Interpreter','latex','FontSize',24);

    xlabel('Position along the string (m)','Interpreter','latex','FontSize',22);

    ylabel('Transverse displacement (m)','Interpreter','latex','FontSize',22);

    legend(string_plot,'String and masses','Interpreter','latex','Location','northeast');

    grid on
    axis equal
    xlim([-0.05*string_length 1.05*string_length])
    ylim([-y_limit y_limit])
    set(gca,'XLimMode','manual','YLimMode','manual')

    video = VideoWriter('vibrating_string.mp4','MPEG-4');
    video.FrameRate = 30;
    video.Quality = 100;
    open(video)

    frame_step = 20;

    for i = 1:frame_step:length(tlist)
        set(string_plot,'YData',Ulist(:,i));
        drawnow
        frame = getframe(fig);
        writeVideo(video,frame);
    end

    close(video)

end