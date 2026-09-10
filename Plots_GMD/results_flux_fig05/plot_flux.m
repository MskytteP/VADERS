clear
close all

load("ADE1d_CASE1and2_3xPE.mat","t_crit_rob","t_crit_dir");

load("ADE1d_var_alpha_100.mat")
Flux_air_100 = Flux_air_cell{1};
Flux_surf_100 = Flux_surf_cell{1};
air_mass_100 = air_mass_cell{1};

load("ADE1d_var_alpha_0.1.mat")
Flux_air_01 = Flux_air_cell{1};
Flux_surf_01 = Flux_surf_cell{1};
air_mass_01 = air_mass_cell{1};

load("ADE1d_var_alpha_0.01.mat")
Flux_air_001 = Flux_air_cell{1};
Flux_surf_001 = Flux_surf_cell{1};
air_mass_001 = air_mass_cell{1};

a_crit = 0.01;
a_crit_arr = repelem(a_crit,length(t_cell{1}));

frac = [0.01 0.1 100];

Flux_air_cell = {};
Flux_air_cell{1} = Flux_air_001;
Flux_air_cell{2} = Flux_air_01;
Flux_air_cell{3} = Flux_air_100;

Flux_surf_cell = {};
Flux_surf_cell{1} = Flux_surf_001;
Flux_surf_cell{2} = Flux_surf_01;
Flux_surf_cell{3} = Flux_surf_100;

air_mass_cell = {};
air_mass_cell{1} = air_mass_001;
air_mass_cell{2} = air_mass_01;
air_mass_cell{3} = air_mass_100;



%% Figure

fig = figure();
set(fig,"position",[1 1 800 480]);

%% Panel (a)

subplot(2,2,1)

handles = [];
labels  = {};

for i = 1:2
    h = plot(t_cell{1},air_mass_cell{i}, ...
        "Linewidth",3);
    hold on

    handles(end+1) = h;
    labels{end+1} = sprintf( ...
        '$\\Lambda_{\\textrm{res}}^*/v_{\\textrm{dep}}^* = %g$', ...
        frac(i));
end

plot(t_cell{1},a_crit_arr,"--k")

yl = ylim();

line([t_crit_rob(2) t_crit_rob(2)],yl, ...
    "color","r","linewidth",2);

line([t_crit_dir(2) t_crit_dir(2)],yl, ...
    "color","b","linewidth",2);

legend(handles,labels, ...
       "interpreter","latex", ...
       "fontsize",16)

ylabel("Airborne mass", ...
       "fontsize",18)

xlabel("$t$", ...
       "interpreter","latex", ...
       "fontsize",20)

title("(a)", ...
      "interpreter","latex", ...
      "fontsize",20)

xlim([0.4 3])
ylim([0 0.1])

set(gca,"fontsize",18)

%% Panel (b)

subplot(2,2,2)

plot(t_cell{1},air_mass_cell{3}, ...
     "linewidth",3, ...
     "color",[0.47 0.67 0.19]);

hold on

plot(t_cell{1},a_crit_arr,"--k")

yl = ylim();

line([t_crit_rob(2) t_crit_rob(2)],yl, ...
    "color","r","linewidth",2);

line([t_crit_dir(2) t_crit_dir(2)],yl, ...
    "color","b","linewidth",2);

legend(sprintf('$\\Lambda_{\\textrm{res}}^*/v_{\\textrm{dep}}^* = %g$', ...
       frac(3)), ...
       "interpreter","latex", ...
       "fontsize",16)

ylabel("Airborne mass", ...
       "fontsize",18)

xlabel("$t$", ...
       "interpreter","latex", ...
       "fontsize",20)

title("(b)", ...
      "interpreter","latex", ...
      "fontsize",20)

xlim([0.4 3])
ylim([0 0.5])

set(gca,"fontsize",18)

%% Panel (c)

subplot(2,2,[3 4])

for i = 1:2

    plot(t_cell{1}, ...
        -(Flux_air_cell{i}-Flux_surf_cell{i}), ...
        "linewidth",3);

    hold on

end

ylabel("$J_{\\textrm{surf}\\rightarrow\\textrm{air}}$", ...
       "interpreter","latex", ...
       "fontsize",18)

xlabel("$t$", ...
       "interpreter","latex", ...
       "fontsize",20)

title("(c)", ...
      "interpreter","latex", ...
      "fontsize",20)

xlim([0.5 2.8])
ylim([-0.1 0.05])

set(gca,"fontsize",18)

