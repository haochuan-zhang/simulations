clear all
clc
close all

dir_name = 'H_Unif0_Discr1_256'
%cd(dir_name) % folder name
fig_list= ls('*.fig') % m * n char, where first two are '.' and '..'
%fig_list = fig_list([3:end],:); %remvoe the '.' and '..'
[fig_num, tmp] = size(fig_list);
for kk= 1:fig_num 
    close all
    fig_list(kk, :)
    end_idx = strfind(fig_list(kk, :), '.fig')+3;
    open(fig_list(kk, 1:end_idx));
    h_line=get(gca,'Children');%get line handles
    xdata=get(h_line,'Xdata');
    ydata=get(h_line,'Ydata');
    tags=get(h_line,'DisplayName')
%     'AMP-SE'
%     'VAMP-SE'
%     'AMP-Algo'
%     'VAMP-Algo'

    AMP_SE(kk,:) = ydata{1,1};
    VAMP_SE(kk,:) = ydata{2,1};
    AMP_Algo(kk,:) = ydata{3,1};
    VAMP_Algo(kk,:) = ydata{4,1};    
end
save([dir_name]);
