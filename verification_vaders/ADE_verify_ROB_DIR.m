function results = ADE_verify_ROB_DIR(D,N)
% Input is D: Turbulent diffusion constant and N: Number of cells in numerical solver
################################
# FUNCTIONS FOR ANALYTICAL SOLUTIONS
################################

function eigensval= newton_rhapson_tan(Pe,low,upp,N,option)
  search_xi = linspace(low,upp,N);  % array for input values. Get sorted later
  err = 1e-8; % err value used to convergence
  dif = 1;  % initial error value before starting loop

  first_eig = zeros(1,length(search_xi));


  for l = 1:length(search_xi)
    dif = 1;  % initial error value before starting loop
    xi = search_xi(l);
    i = 1;  % used for counting if the expression does not converge
    while dif > err && i < 10000
      xi_prev = xi;

      if option == "pelte2"

        fxi = Pe/2*sin(xi) - xi*cos(xi);
        dfxi = cos(xi)*(Pe/2-1) + xi*sin(xi);

      elseif option == "x"

        fxi = xi*cosh(xi) - Pe/2*sinh(xi);
        dfxi = cosh(xi)*(1-Pe/2) + xi*sinh(xi);

      else
        printf("choose either lte2 or gt2")
      end

      xi = xi - fxi/dfxi;

      dif = abs(xi_prev - xi);
      i+=1;
      if i == 9999
        xi = 0
        break
      endif
    end
    first_eig(l) = xi;

  end
  %% Sorting values to get uniques positive

  eigen_vals = unique(sort(first_eig));
  tjek = eigen_vals(2:end)-eigen_vals(1:end-1) < 1e-4;
  vals = eigen_vals(tjek);
  pos = find(eigen_vals>1e-10);
  eigen_vals = eigen_vals(pos);
  for jj = 1:length(eigen_vals(tjek))
    kk = find(eigen_vals == vals(jj));  %real
    eigen_vals(kk) = 0;
  endfor

  eigen_vals = eigen_vals(find(0<eigen_vals));  %real
  eigensval = eigen_vals(find(eigen_vals<upp));


end

function areturn = ADE_a(z, Pe, xi, a_n, t)
  v = zeros(1,length(z));

  if Pe <=2

    for nmode = 1:length(xi)
      v = v + a_n(nmode).*(cos(xi(nmode).*z) - Pe/(2*xi(nmode))*sin(xi(nmode)*z)).*exp(-1/Pe*xi(nmode).^2*t);
    endfor

    areturn = v.*exp(-Pe/2*z-Pe/4*t);

   elseif Pe > 2

    v = a_n(1).*(cosh(xi(1).*z) - Pe/(2*xi(1))*sinh(xi(1).*z)).*exp(1/Pe*xi(1).^2*t);
    for nmode = 2:length(xi)
      v = v + a_n(nmode).*(cos(xi(nmode).*z) - Pe/(2*xi(nmode))*sin(xi(nmode)*z)).*exp(-1/Pe*xi(nmode).^2*t);
    endfor

    areturn = v.*exp(-Pe/2*z-Pe/4*t);
  endif

endfunction

function integ=num_integral(xi,z, Pe, sigma, A)


   fun1 = @(z) (cos(xi.*z) - Pe./(2.*xi).*sin(xi.*z)).*exp(Pe/(2).*z).*A.*exp(-(z-1/2).^2/(2*sigma^2));
   integ_taller = integral(fun1,0,1);
   fun2 = @(z) (cos(xi.*z) - Pe/(2.*xi)*sin(xi.*z)).^2;
   integ_navner = integral(fun2,0,1);

   integ = integ_taller/integ_navner;
end

function integ_hyp = num_integral_hyp(xi, z, Pe, sigma, A)

  fun1 = @(z) (cosh(xi.*z) - Pe./(2.*xi).*sinh(xi.*z)).*exp(Pe/(2).*z).*A.*exp(-(z-1/2).^2/(2*sigma^2));
  integ_taller = integral(fun1,0,1);
  fun2 = @(z) (cosh(xi.*z) - Pe/(2.*xi)*sinh(xi.*z)).^2;
  integ_navner = integral(fun2,0,1);

  integ_hyp = integ_taller/integ_navner;
endfunction

%% Solver parameters %%

save_n = 10000;       % save every save_nth step


%% Parameters %%

L = 100;    % Domain length [m]

v_set = 3e-3  % Settling velocity due to gravity[m/s]  % Reference aerobiology book
pe = v_set*L/D  % Peclet number




v_dep = 0   % Surface deposition [m/s]
lambda_r = 0;  % Resuspension rate [1/s]
alpha = 0; % Inactivation rate constant [1/s]

%% Discretization %%
dz = L/N;   % Cell size
z = linspace(dz/2,L-dz/2,N);

z = z/L;        % nondimensionalization
dz = z(2)-z(1);   % nondimensionalization

CFL_adv  = 0.5;              % advection CFL <= 1
CFL_diff = 0.25;             % diffusion CFL <= 0.5
dt1 = CFL_adv * dz / max(abs(v_set), 1e-16);
dt2 = CFL_diff * dz^2 / D;
dt = min(dt1, dt2)*0.1;


%% initial conditions %%
% --- AIR ---
sigma =  0.05;  % spread
A = 1/(sigma*sqrt(2*pi)); % amplitude
a_ini = A*exp(-(z-0.5).^2/(2*sigma^2)); % initial distribution
a = a_ini';

% --- SURFACE ---
a_s = 0;
b_s = 0;    % Removal of material

%% Accounting for mass

mass_air = sum(a)*dz
mass_ini = mass_air + a_s + b_s
mass_fraction_stop = 0.01;

%% For plotting %%
t_arr = [];
as_arr = [];
bs_arr = [];
Mass_arr = [];
L2rel_arr = [];
L2_norm_arr = [];
L2_err_arr = [];

o=0;

add_surf = v_dep*L^2/(dz*L*D);
resuspension = lambda_r*L^2/D;
removal = alpha*L^2/D;


### SETTING UP ANALYTICAL #########

ini_arr_ana = zeros(length(z),1);

fail_beta = [];
fail_index = [];
ini_fail = [];
u = 1;

  if pe <=2
    xi = newton_rhapson_tan(pe,0,220,1000,"pelte2");
    pe
  elseif pe >2
    xi_1 = newton_rhapson_tan(pe,0,50000,4000,"x");%0.3
    xi_rest = newton_rhapson_tan(pe,0,220,1000,"pelte2");
    xi = [xi_1 xi_rest];

  end

  if pe < 2 && xi(1) < 0.1
    xi = xi(2:end);
  end

  a_n = zeros(1,length(xi));

   if pe <=2
     for k = 1:length(xi)
        a_n(k) = num_integral(xi(k),z, pe, sigma, A); % forced beta
      end

    elseif pe >2
      a_n(1) = num_integral_hyp(xi(1),z,pe,sigma,A);
      for k = 2:length(xi)
        a_n(k) = num_integral(xi(k),z, pe, sigma, A);
      endfor
    endif

  if length(xi) == 70%126
    a_n_arr = a_n;
    xi_arr = xi;
  elseif length(xi) == 71%127;
    a_n_arr = a_n(1:end-1);
    xi_arr = xi(1:end-1);
  end
  if trapz(z,ADE_a(z,pe,xi,a_n,0)) > 0.005
    ini_arr_ana = ADE_a(z,pe,xi,a_n,0);
  else
    ini_arr_ana(:) = 0;
    fail_beta(u) = v_set;
    ini_fail(:,u) = ADE_a(z,pe,xi,a_n,0);
    u+=1;
  end




n = 0;
middle_saved = false;

while true

  n += 1;
  %% Fluxes computed at interfaces (N+1 faces)
  F = zeros(N+1,1);



  %% Bottom Boundary Flux
  J_surf = add_surf*a(1) - resuspension*a_s;
  F(1) = -J_surf;

  %% internal fluxes
  F(2:N) = -a(2:N) - (1/pe)*(a(2:N)-a(1:N-1))/dz;


  %% Top Boundary Flux
 % F(N+1) = 0;    % Closed system
  a_top = 0; % Dirichlet

  a_up = 0;

  F_adv_top = - a_up;
  F_diff_top = -1/pe*(a_top-a(N)) / dz;

  F(N+1) = F_adv_top + F_diff_top;

  %% Update air concentration
  a = a + (dt/dz)*(F(1:N) - F(2:N+1));

  %% Update surface concentration
  removal_surf = removal*a_s;
  a_s = a_s + (J_surf - removal_surf)*dt;
  b_s = b_s + removal_surf*dt;

  %% Track mass
  mass_air = sum(a)*dz;
  total_mass = mass_air + a_s + b_s;

  mass_fraction = mass_air/mass_ini;

  if !middle_saved && mass_fraction < 0.5

    a_num_middle = a;
    a_ana_middle = ADE_a(z,pe,xi_arr,a_n_arr,n*dt);

    tau_middle = n*dt;

    middle_saved = true;

  endif

  if mass_fraction < mass_fraction_stop

      fprintf('\n');
      fprintf('Stopping criterion reached.\n');
      fprintf('Remaining airborne mass = %.4e\n',mass_fraction);
      tau_final = n*dt;
      mass_fraction_final = mass_fraction;
      break

  endif

  if rem(n,save_n) == 0
    o+=1;
    t_arr(o) = n*dt;
    as_arr(o) = a_s;
    Mass_arr(o) = total_mass;



    a_ana = ADE_a(z,pe,xi_arr,a_n_arr,t_arr(o));

    L2err = sqrt(sum((a(:)-a_ana(:)).^2)*dz);

    L2rel = L2err / sqrt(sum(a_ana(:).^2)*dz);
    L2rel_arr(o) = L2rel;

    L2_err_arr(o) = L2err;
    L2_norm_arr(o) = sqrt(sum(a_ana(:).^2)*dz);

    mass_ana = sum(a_ana)*dz;

    rel_dif = (mass_ana - Mass_arr(o))/mass_ana*100;
    if rem(n,save_n*10) == 0
      fprintf('pe = %d | tau = %3.2f | remaining mass = %6.2f %%\n', pe, n*dt, 100*mass_fraction);
    endif

%{
    subplot(2,2,[1 2])
    plot(a,z, "Linewidth",2, '*r')
    hold on
    plot(a_ana,z, "Linewidth",2, 'ok')
    xlabel("a")
    ylabel("z")
    xlim([0 10])
    legend("Numerical", "Analytical")
    titl = sprintf("Concentration profile at tau = %02d (t = %02d hours)",t_arr(o),t_arr(o)*L/v_set/3600);
    title(titl)
    hold off
    set(gca,"FontSize", 20)

    subplot(2,2,3)
    hold on
    plot(t_arr(o),L2rel_arr(o),'*r',"Markersize",10)
    ylabel("Relative L_2 error")
    xlabel("\\tau")
    %xlim([0 t_end])
    hold off
    legend("Relative difference from analytical")
    set(gca,"FontSize", 20)

    subplot(2,2,4)
    hold on
    plot(t_arr(o),L2_err_arr(o),'*r',"Markersize",10)
    plot(t_arr(o), L2_norm_arr(o), '+b',"Markersize",10)
    ylabel("Relative L_2 error")
    legend("L2 error","||a_{ana}||_2")
    xlabel("\\tau")
    %xlim([0 t_end])
    hold off
    set(gca,"FontSize", 20)

    pause(0.001)
%}
  endif
end

results = struct();

% Case information
results.Pe = pe;
results.N = N;
results.D = D;
results.L = L;
results.dt = dt;
results.tau_final = tau_final;
results.mass_fraction_final = mass_fraction_final;

% Numerical settings
results.CFL_adv = CFL_adv;
results.CFL_diff = CFL_diff;

% Grid
results.z = z;
results.dz = dz;

% Initial conditions
results.a_initial = a_ini;
results.a_ana_initial = ADE_a(z,pe,xi_arr,a_n_arr,0);

% Time series
results.tau = t_arr;
results.L2rel = L2rel_arr;
results.L2err = L2_err_arr;
results.L2norm = L2_norm_arr;

results.mass_num = Mass_arr;

% Mass 0.5 state
results.a_num_middle = a_num_middle;
results.a_ana_middle = a_ana_middle;
results.tau_middle = tau_middle;

% Final state
results.a_num_final = a;
results.a_ana_final = a_ana;
results.a_s_final = a_s;
results.b_s_final = b_s;

a_ana_final = ADE_a(z,pe,xi_arr,a_n_arr,tau_final);
L2err_final = sqrt(sum((a(:)-a_ana_final(:)).^2)*dz);
L2norm_final = sqrt(sum(a_ana_final(:).^2)*dz);
L2rel_final = L2err_final / L2norm_final;

results.L2rel_final = L2rel_final;
results.L2err_final = L2err_final;
results.L2norm_final = L2norm_final;






% Analytical expansion
results.xi = xi_arr;
results.a_n = a_n_arr;

end

%filename = sprintf('verification_Pe_%0.2f_N_%d_tend_%g.mat', pe, N, t_end);
%save(filename,'results')



