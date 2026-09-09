function savefil = vaders_main(frac)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates the figure with different air-clearance times
% for a specific (frac = lambda^*_res/v^*_dep) and constant
% diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

D = [0.075 0.3 1];   % Turbulent Diffusion coefficient (Kz)
L = 100;    % Domain length [m]

v_set = 3e-3;%3.01e-5  % Settling velocity due to gravity[m/s]

Pe_arr = v_set*L./D;  % Peclet numbers

v_dep = 5e-4;  % Surface deposition velocity [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Discretization %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

N = 250;    % number of cells
dz = L/N;   % Cell size
z = linspace(dz/2,L-dz/2,N);

z = z/L;        % nondimensionalization
dz = z(2)-z(1);   % nondimensionalization

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Surface dynamics parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

v_dep_star = v_dep*L/(dz*L*v_set);
lambda_r = v_dep_star*frac/(dz*L); % it is lambda^*_r: frac is the argument in function (Lambda^*_res/v^*_dep)
frac_show = lambda_r/v_dep_star*dz*L

alpha_arr = [0 5e-6 5e-5 5e-4 5e-3]; % inactivation rate constants [1/s]

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% initial conditions %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% AIR  %%%%%
sigma =  0.05;  % spread
A = 1/(sigma*sqrt(2*pi)); % amplitude
a_ini = A*exp(-(z-0.5).^2/(2*sigma^2)); % initial distribution
a = a_ini';

%%%% SURFACE %%%%
a_s=0;    % Surface Concentration
b_s=0;    % Removal from surface (for accounting)

%%%%% Accounting for mass %%%%%%%

mass_air = sum(a)*dz
mass_ini = mass_air + a_s + b_s
a_crit = 0.01*mass_ini;
total_mass = mass_ini;

%%%%%%%% For plotting %%%%%%%%%%


Flux_top_cell = {};
Flux_air_cell = {};
Flux_surf_cell = {};
as_cell = {};
mass_cell = {};


t_crit_num = zeros(length(Pe_arr),length(alpha_arr));

for h = 1:length(alpha_arr)

  alpha = alpha_arr(h)

  for j = 1:length(Pe_arr)
    dt_max = 0.5 * (dz)^2 / (D(j) + 0.5*v_set*dz);
    gamma = 0.2;
    dt = gamma * dt_max
    b_s = 0;    % accounting removal of material from inactivation
    add_surf = v_dep_star
    resuspension = lambda_r
    removal = alpha*L/v_set;    % alpha^*

    pe = Pe_arr(j)
    total_mass = mass_ini
    t_arr = []; % zeros(Nt,1);
    as_arr = [];  % zeros(Nt,1);
    bs_arr = [];
    Mass_arr = [];  % zeros(Nt,1);

    Flux_top_arr = [];
    Flux_air_arr = [];
    Flux_surf_arr = [];

    a = a_ini';
    a_s = 0;
    o=0;
    k=1;
    while total_mass > a_crit
      %% Fluxes computed at interfaces (N+1 faces)
      F = zeros(N+1,1);



      %% Bottom Boundary Flux
      J_surf = add_surf*a(1) - resuspension*a_s;
      F(1) = -J_surf;
      Flux_air = add_surf*a(1);
      Flux_surf = resuspension*a_s;

      %% internal fluxes
      F(2:N) = -a(2:N) - (1/pe)*(a(2:N)-a(1:N-1))/dz;

      %% Top Boundary Flux

      % Closed system %
      % F(N+1) = 0;

      % Dirichlet %
      a_top = 0;

      a_up = 0;

      F_adv_top = - a_up;
      F_diff_top = -1/pe*(a_top-a(N)) / dz;

      F(N+1) = F_adv_top + F_diff_top;
      Flux_top = F(N+1);
      %% Update air concentration
      a = a + (dt/dz)*(F(1:N) - F(2:N+1));

      %% Update surface concentration
      removal_surf = removal*a_s;
      a_s = a_s + (J_surf - removal_surf)*dt;
      b_s = b_s + removal_surf*dt;

      %% Track mass
      mass_air = sum(a)*dz;
      total_mass = mass_air;

      if rem(k,5000) == 0
        o+=1;
        t_arr(o) = k*dt;
        as_arr(o) = a_s;
        Mass_arr(o) = total_mass;
        Flux_top_arr(o) = Flux_top;
        Flux_air_arr(o) = Flux_air;
        Flux_surf_arr(o) = Flux_surf;


%%%%%%%% PLOT %%%%%%%%%%%
        subplot(2,2,[1 2])
        plot(a,z, "Linewidth",3, 'b')
        xlabel("a")
        ylabel("z")
        xlim([0 10])
        legend("Numerical")
        titl = sprintf("Concentration profile at tau = %02d (t = %02d hours)",t_arr(o),t_arr(o)*L/v_set/3600);
        title(titl)
        hold off
        set(gca,"FontSize", 20)


        subplot(2,2,3)
        hold on
        plot(t_arr(o),a_s,'*r','Markersize', 10)
        xlabel("\\tau")
        ylabel("a_s")
        xlim([0 20])
        ylim([0 mass_ini])
        legend("Surface concentration")
        hold off
        set(gca,"FontSize", 20)

        subplot(2,2,4)
        hold on
        plot(t_arr(o),b_s,'*r','Markersize', 10)
        xlabel("\\tau")
        ylabel("b_s")
        xlim([0 20])
        ylim([0 mass_ini])
        legend("Surface removal")
        hold off
        set(gca,"FontSize", 20)
        pause(0.01)

    endif

      k+=1;

    end
    t_crit_num(j,h) = k*dt
    Flux_top_cell{j,h} = Flux_top_arr;
    Flux_air_cell{j,h} = Flux_air_arr;
    Flux_surf_cell{j,h} = Flux_surf_arr;
    as_cell{j,h} = as_arr;
    mass_cell{j,h} = Mass_arr;

    t_cell{j,h} = t_arr;

end

end
    savefil = sprintf("results/ADE1d_var_alpha_%d.mat",frac);
    save(savefil, "Pe_arr", "t_crit_num","v_set","L","D","lambda_r","dz","v_dep","v_dep_star","alpha_arr","t_cell","Flux_top_cell", "Flux_air_cell","Flux_surf_cell")
end




