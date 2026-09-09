
function results = ADE_verify_masscons(frac)
%% Parameters %%

save_n = 10000;

D = 0.075% [0.075 0.3 1];%0.08 0.1 0.2  0.06linspace(0.05,0.5,10);%0.12105;    % Diffusion coefficient
L = 100;    % Domain length [m]

v_set = 3e-3;%3.01e-5  % Settling velocity due to gravity[m/s]  % Reference aerobiology book

Pe_arr = v_set*L./D;
%alpha = D;    % alpha is diffusion in old solution


%lambda_arr = [5e-8 8e-8];
v_dep = 5e-4; %ARG5e-3;%0.0001;   % Surface deposition [m/s] Reference Cheolwoon 2018 fungi measurements


%% Discretization %%
N = 250;    % number of cells
dz = L/N;   % Cell size
z = linspace(dz/2,L-dz/2,N);

z = z/L;        % nondimensionalization
dz = z(2)-z(1);   % nondimensionalization


v_dep_star = v_dep*L/(dz*L*v_set);
lambda_r = v_dep_star*frac/(dz*L);%
frac_show = lambda_r/v_dep_star*dz*L

alpha_arr = [0];


%% initial conditions %%
% --- AIR ---
sigma =  0.05;  % spread
A = 1/(sigma*sqrt(2*pi)); % amplitude
a_ini = A*exp(-(z-0.5).^2/(2*sigma^2)); % initial distribution
a = a_ini';

% --- SURFACE ---
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

%o=0;



t_crit_num = zeros(length(Pe_arr),length(alpha_arr));

for h = 1:length(alpha_arr)

  alpha = alpha_arr(h)

  for j = 1:length(Pe_arr)
    dt_max = 0.5 * (dz)^2 / (D(j) + 0.5*v_set*dz);
    gamma = 0.2;
    dt = gamma * dt_max
    b_s = 0;    % Removal of material
    add_surf = v_dep_star
    resuspension = lambda_r
    removal = alpha*L/v_set;

    pe = Pe_arr(j)
    total_mass = mass_ini
    t_arr = [];%zeros(Nt,1);
    mass_air_arr = [];
    mass_surface_arr = [];
    mass_removed_arr = [];
    mass_total_arr = [];

    a = a_ini';
    a_s = 0;
    o=0;
    k=1;

    tau_final = 50;
    while true
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
      F(N+1) = 0;    % Closed system
      %{
      a_top = 0;

      a_up = 0;%a(N);

      F_adv_top = - a_up;
      F_diff_top = -1/pe*(a_top-a(N)) / dz;

      F(N+1) = F_adv_top + F_diff_top;
      Flux_top = F(N+1);
      %}

      %% Update air concentration
      a = a + (dt/dz)*(F(1:N) - F(2:N+1));

      %% Update surface concentration
      %J_surf = add_surf*a(1) - resuspension*a_s;
      removal_surf = removal*a_s;
      a_s = a_s + (J_surf - removal_surf)*dt;
      b_s = b_s + removal_surf*dt;

      %% Track mass
      mass_air = sum(a)*dz;


      if rem(k,save_n) == 0
        o+=1;
        t_arr(o) = k*dt;

        mass_air_arr(o) = mass_air;
        mass_surface_arr(o) = a_s;
        mass_removed_arr(o) = b_s;
        mass_total_arr(o) = mass_air + a_s + b_s;



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
        plot(t_arr(o),mass_total_arr(o)-mass_ini,'*r','Markersize', 10)
        xlabel("\\tau")
        ylabel("Airborne mass")
        %xlim([0 20])
        %ylim([0 mass_ini])
        legend("Surface removal")
        hold off
        set(gca,"FontSize", 20)
        pause(0.01)

      endif

    tau = k*dt;
    if tau >=tau_final || mass_air<=a_crit
      tau_final = tau;
      break
    endif

      k+=1;

    end
results = struct();

%% Case setup
results.D = D(j);
results.Pe = pe;

results.L = L;
results.N = N;

results.v_set = v_set;
results.v_dep = v_dep;

results.lambda_r = lambda_r;
results.alpha = alpha;

results.frac = frac;

%% Numerical setup
results.dz = dz;
results.dt = dt;

%% Initial condition
results.mass_ini = mass_ini;
results.a_crit = a_crit;

%% Time series
results.tau = t_arr;

results.mass_air = mass_air_arr;
results.mass_surface = mass_surface_arr;
results.mass_removed = mass_removed_arr;
results.mass_total = mass_total_arr;

%% Final values
results.mass_air_final = mass_air;
results.mass_surface_final = a_s;
results.mass_removed_final = b_s;
results.mass_total_final = mass_air + a_s + b_s;

%% Conservation error
results.mass_error = mass_total_arr - mass_ini;
results.relative_mass_error = (mass_total_arr - mass_ini)/mass_ini;

%% Stopping time
results.tau_final = k*dt;


end

end
end




