function time_hist=choose_time_hist_vs_single
% Construct a questdlg with three options
choice = questdlg('This file contains a single solution or a time history?',...
    'Type of output file',...
    'Single solution','Time history','Single solution');
% Handle response
switch choice
    case 'Single solution'
        disp([choice ' ok: single solution'])
        time_hist = 0;
    case 'Time history'
        disp([choice ' ok: time history'])
        time_hist = 1;
end
end