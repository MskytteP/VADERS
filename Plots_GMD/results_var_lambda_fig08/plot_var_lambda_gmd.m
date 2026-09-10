clear
close all

files = dir("ADE1d_var_t_lambda_*.mat");
file_names = {files.name};

hour = zeros(length(file_names),1);

for i = 1:length(file_names)
  fname = files(i).name;
  token= regexp(fname,'ADE1d_var_t_lambda_(.*)\.mat','tokens');
  hour(i) = str2double(token{1}{1});

  load(file_names{i});
  t_name = sprintf("t_crit_%d",i);
  time_name = sprintf("t_cell_%d",i);
  flux_s_name = sprintf("flux_s_%d",i);
  eval([t_name " =  t_crit_num;"]);
  eval([time_name " = t_cell;"]);
  eval([flux_s_name " = Flux_surf_cell;"]);
end


[hour,idx] = sort(hour);

hour_tau = hour*v_set/L*3600;

for i = 1:length(file_names)
    varname = sprintf("t_crit_%d", i);
    temp = eval(varname);
    arr_pe4(1,i) = temp(1,:);
end

for i = 1:length(file_names)
    varname = sprintf("t_crit_%d", i);
    temp = eval(varname);
    arr_pe1(1,i) = temp(2,:);
end

for i = 1:length(file_names)
    varname = sprintf("t_crit_%d", i);
    temp = eval(varname);
    arr_pe03(1,i) = temp(3,:);
end

arr_pe03 = arr_pe03(idx);
arr_pe1 = arr_pe1(idx);
arr_pe4 = arr_pe4(idx);

figure("position",[1 1 800 480])
%setNicePlotDefaults()
set(gca,'ColorOrderIndex',1)

plot(hour_tau(1:2:end),arr_pe4(1:2:end),'*-','Linewidth',2)
hold on
plot(hour_tau(1:2:end),arr_pe1(1:2:end),'^-','Linewidth',2)
plot(hour_tau(1:2:end),arr_pe03(1:2:end),'d-','Linewidth',2)

legend('$Pe=4$','$Pe=1$','$Pe=0.3$','interpreter','latex', "fontsize",18)
set(legend,'location', "northwest")

xlabel('$t_0$','interpreter','latex', "fontsize",30)
ylabel('$t_c$','interpreter','latex', "fontsize",30)
set(gca,"Fontsize",20)
%tit=sprintf('Relationship of $t_c$ and $\\Lambda_{\\text{res}}(t)$');
%set(tit, "fontsize",25)

%title(tit,'interpreter','latex')

