clear
close all

files = dir("ADE1d_var_diff_*.mat");
file_names = {files.name};

frac = zeros(length(file_names),1);

for i = 1:length(file_names)
  fname = files(i).name;
  token= regexp(fname,...
      'ADE1d_var_diff_(.*)\.mat',...
      'tokens');

  frac(i) = str2double(token{1}{1});

  load(file_names{i});
  t_name = sprintf("t_crit_%d",i);
  time_name = sprintf("t_cellpe4_%d",i);
  eval([t_name " =  t_crit_num;"]);
  eval([time_name " = t_cell;"]);
end

arr_pe4 = zeros(length(alpha_arr),length(file_names));
t_pe4 = cell(length(alpha_arr), length(file_names));

for j = 1:length(alpha_arr)
  for i = 1:length(file_names)
      varname = sprintf("t_crit_%d", i);
      tvarname = sprintf("t_cellpe4_%d",i);
      temp = eval(varname);
      tempt = eval(tvarname);
      arr_pe4(:,i) = temp(1,:);
      t_pe4{j,i} = tempt{1,j};
  end
end

[frac,idx] = sort(frac);
arr_pe4 = arr_pe4(:,idx);
t_pe4 = t_pe4(:,idx);



set(gca,'ColorOrderIndex',1)
figure("position",[1 1 800 480])

hpe4=loglog(frac,arr_pe4,'*-',"Linewidth",3,"Markersize",8,"DisplayName","Pe4")


legend(hpe4,...
{'$\alpha^* = 0$',...
 '$\alpha^* = 1.67\times10^{-1}$',...
 '$\alpha^* = 1.67$',...
 '$\alpha^* = 1.67\times10^{1}$',...
 '$\alpha^* = 1.67\times10^{2}$'},...
'location','northwest',...
'interpreter','latex')

xlabel('$\Lambda^*_{\text{res}}/v^*_{\text{dep}}$','interpreter','latex')
ylabel('$t_{\text{c}}$','interpreter','latex')
title("(b)",'interpreter','latex')
set(gca,"Fontsize",24)
