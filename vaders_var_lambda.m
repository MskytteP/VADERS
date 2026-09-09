
function savefil = vaders_var_lambda(hours)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates the figure with different air-clearance times
% for a specific number of hours with high resuspension and with
% varying lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters %%

D = [0.075 0.3 1]; % Turbulent diffusion coefficient
L = 100;   % Domain length [m]

v_set = 3e-3; % Settling velocity due to gravity[m/s]
Pe_arr = v_set*L./D;  % Peclet numbers
v_dep = 5e-4 % Surface deposition [m/s] Reference Cheolwoon 2018 fungi measurements


%% Discretization %%
N = 250;    % number of cells
dz = L/N;   % Cell size
z = linspace(dz/2,L-dz/2,N);

z = z/L;        % nondimensionalization
dz = z(2)-z(1);   % nondimensionalization

lambda_min = 5e-8;  % [1/s]
lambda_max = 1e-4; % Heavy vehicle traffic(~5e-5) [1/s]
alpha_arr = [0];  % inactivation rate constant [1/s]



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
lambda_cell = {};

t_crit_num = zeros(length(Pe_arr),length(alpha_arr));

for h = 1:length(alpha_arr)

  alpha = alpha_arr(h)

  for j = 1:length(Pe_arr)
    dt_max = 0.5 * (dz)^2 / (D(j) + 0.5*v_set*dz);
    gamma = 0.2;
    dt = gamma * dt_max
    b_s = 0;    % Accounting for removal of material
    add_surf = v_dep*L/(dz*L*v_set) % v^*_dep
    resuspension_max = lambda_max*L/v_set % dimless
    resuspension_min = lambda_min*L/v_set % dimless
    removal = alpha*L/v_set;  % alpha^*

    pe = Pe_arr(j)
    total_mass = mass_ini
    t_arr = [];
    as_arr = [];
    bs_arr = [];
    Mass_arr = [];

    Flux_top_arr = [];
    Flux_air_arr = [];
    Flux_surf_arr = [];
    resuspension_arr = [];

    tau_0 = v_set/L*3600*hours; % hours of traffic

    a = a_ini';
    a_s = 0;
    o=0;
    k=1;
    while total_mass > a_crit
      %% Fluxes computed at interfaces (N+1 faces)
      F = zeros(N+1,1);

      resuspension = lambda_min*L/v_set + (lambda_max*L/v_set-lambda_min*L/v_set)./(1+exp(30*(k*dt-tau_0)));


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
        resuspension_arr(o) = resuspension;

        %%%% PLOTTING %%%%%

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
        xlim([0 10])
        ylim([0 mass_ini])
        legend("Surface concentration")
        hold off
        set(gca,"FontSize", 20)

        subplot(2,2,4)
        hold on
        plot(t_arr(o),resuspension_arr(o),"*r", "Markersize",10)
        xlabel("\\tau")
        ylabel("\\Lambda_{res}")
        xlim([0 10])
        ylim([0 resuspension_max])
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
    lambda_cell{j,h} = resuspension_arr;
    t_cell{j,h} = t_arr;

end

end
    savefil = sprintf("results/ADE1d_var_t_lambda_%d.mat",hours);
    save(savefil, "tau_0", "Pe_arr", "t_crit_num","v_set","L","D","lambda_max","lambda_min","dz","v_dep","alpha_arr","t_cell","Flux_top_cell", "Flux_air_cell","Flux_surf_cell","lambda_cell")
end




