
function savefil = vaders_var_diff(frac)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates the figure with different air-clearance times
% for a specific (frac = lambda^*_res/v^*_dep) and with varying
% diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters %%
L = 100;    % Domain length [m]

v_set = 3e-3;%3.01e-5  % Settling velocity due to gravity[m/s]
kappa = 0.4; % von karman coefficient
u_fric = 0.3 % friction velocity Open flat terrain, wind 5 m/s

Pe0 = v_set/(kappa*u_fric)

v_dep = 5e-4; % Surface deposition [m/s]

%% Discretization and more parameters %%
N = 250;    % number of cells
dz = L/N;   % Cell size
z = linspace(dz/2,L-dz/2,N);

Pe_z = Pe0*L./z
Kz = kappa*u_fric*z;

z = z/L;        % nondimensionalization
dz = z(2)-z(1);   % nondimensionalization

v_dep_star = v_dep*L/(dz*L*v_set);
lambda_star = v_dep_star*frac; % it is lambda^*_r: frac is the argument in function (Lambda^*_res/v^*_dep)
lambda_res = lambda_star*v_set/L; % Resuspension rate [1/s]

alpha_arr = [0 5e-6 5e-5 5e-4 5e-3];  % inactivation rate constants [1/s]


%% initial conditions %%
% AIR %
sigma =  0.05;  % spread
A = 1/(sigma*sqrt(2*pi)); % amplitude
a_ini = A*exp(-(z-0.5).^2/(2*sigma^2)); % initial distribution
a = a_ini';

% SURFACE %
a_s=0;
b_s=0;

%% Accounting for mass

mass_air = sum(a)*dz
mass_ini = mass_air + a_s + b_s
a_crit = 0.01*mass_ini;
total_mass = mass_ini;

%% For plotting %%


Flux_top_cell = {};
Flux_air_cell = {};
Flux_surf_cell = {};
as_cell = {};
mass_cell = {};


t_crit_num = zeros(length(Pe0),length(alpha_arr));

for h = 1:length(alpha_arr)

  alpha = alpha_arr(h)

  for j = 1:length(Pe0)
    dt_max = (dz)^2 / (max(1./Pe_z));
    gamma = 0.3;
    dt = gamma * dt_max
    b_s = 0;    % accounting for removal of material
    add_surf = v_dep_star % v^*_dep
    resuspension = lambda_star % lambda^*_res
    removal = alpha*L/v_set;  % alpha^*

    pe = Pe0(j)
    total_mass = mass_ini
    t_arr = [];
    as_arr = [];
    bs_arr = [];
    Mass_arr = [];

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
      F_diff_top = -1/Pe0*(a_top-a(N)) / dz;

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

%%% PLOTTING %%%%

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
        xlim([0 1])
        ylim([0 mass_ini])
        legend("Surface concentration")
        hold off
        set(gca,"FontSize", 20)

        subplot(2,2,4)
        hold on
        plot(t_arr(o),(Mass_arr(o)-mass_ini)*100,'*r','Markersize', 10)
        xlabel("\\tau")
        ylabel("total mass")
        xlim([0 1])
        legend("Total mass")
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
    savefil = sprintf("results/ADE1d_var_diff_%g.mat",frac);
    save(savefil, "Pe_z","Pe0","Kz", "t_crit_num","v_set","L","lambda_res","dz","v_dep","alpha_arr","t_cell","Flux_top_cell", "Flux_air_cell","Flux_surf_cell","as_cell","mass_cell")
end




