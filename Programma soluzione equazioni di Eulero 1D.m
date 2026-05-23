%% import & format data
function import_Euler_1D
% import raw data
%[filename,pathname] = uigetfile('results_num_*.dat','Select the MATLAB code file');
[filename,pathname] = uigetfile('num_sol.dat','Select the MATLAB code file');
%basename = extractBetween(filename,'results_num_','.dat');

%analyt_name = "results_ana_" + basename + ".dat";
filename = strcat(pathname,filename);
analyt_name = strcat(pathname,"theo_sol.dat");

time_hist=choose_time_hist_vs_single;

if (time_hist==0)
    t=CLOSE1;
    Euler_1D_raw = importdata(filename);
    analyt_raw = importdata(analyt_name);
    variables_names=transpose(extractBetween(Euler_1D_raw.textdata{2,1},'"','"'));
    i = Euler_1D_raw.textdata(4,2);
    nc = str2double(i);
    % build a table from the imported data
    CFD_data{1} = array2table(Euler_1D_raw.data,'VariableNames',variables_names);
else
    fileID = fopen(filename);
    t=0;
    while ~feof(fileID)
        line = fgetl(fileID);
        a = strncmp(' ZONE T="t=',line,11);
        b = strncmp(' VARIABLES=',line,11);
        if b==1
            variables_names=transpose(extractBetween(line,'"','"'));
        end
        if a==1
            t=t+1;
            timechr(t)=extractBetween(line,'"t=','"');
            time(t) = str2double(extractBetween(line,'"t=','"'));
            nc = str2double(extractBetween(line,'I=','F'));
            C = textscan(fileID,'%f %f %f %f %f %f %f');
            data(:,:,t)=cell2mat(C);
            CFD_data{t} = array2table(data(:,:,t),'VariableNames',variables_names);
        end
    end

    fileID = fopen(analyt_name);
    t=0;
    while ~feof(fileID)
        line = fgetl(fileID);
        a = strncmp(' ZONE T="t=',line,11);
        b = strncmp(' VARIABLES=',line,11);
        if b==1
            variables_names=transpose(extractBetween(line,'"','"'));
        end
        if a==1
            t=t+1;
            timechr(t)=extractBetween(line,'"t=','"');
            time(t) = str2double(extractBetween(line,'"t=','"'));
            nc = str2double(extractBetween(line,'I=','F'));
            C = textscan(fileID,'%f %f %f %f %f %f %f');
            th_data(:,:,t)=cell2mat(C);
            theo_data{t} = array2table(th_data(:,:,t),'VariableNames',variables_names);
        end
    end
    
    tmax=t;
    timemin=min(time);
    timemax=max(time);
    for i=1:tmax-2
        dt(i)=time(i+1)-time(i);
    end
    dtmin = min(dt);
    dtmax = max(dt);
end
if (time_hist==1)
    [t,v] = listdlg('PromptString','Select a time:',...
        'SelectionMode','single',...
        'ListString',timechr)
end
for k=1:tmax
    pmindum(k) = min(CFD_data{1,k}.p);
    pmaxdum(k) = max(CFD_data{1,k}.p);
    rmindum(k) = min(CFD_data{1,k}.rho);
    rmaxdum(k) = max(CFD_data{1,k}.rho);
    Tmindum(k) = min(CFD_data{1,k}.T);
    Tmaxdum(k) = max(CFD_data{1,k}.T);
    umindum(k) = min(CFD_data{1,k}.u);
    umaxdum(k) = max(CFD_data{1,k}.u);
    Mmindum(k) = min(CFD_data{1,k}.M);
    Mmaxdum(k) = max(CFD_data{1,k}.M);
    Smindum(k) = min(CFD_data{1,k}.S);
    Smaxdum(k) = max(CFD_data{1,k}.S);
end
pmin = min(pmindum);
pmax = max(pmaxdum);
rmin = min(rmindum);
rmax = max(rmaxdum);
Tmin = min(Tmindum);
Tmax = max(Tmaxdum);
umin = min(umindum);
umax = max(umaxdum);
Mmin = min(Mmindum);
Mmax = max(Mmaxdum);
Smin = min(Smindum);
Smax = max(Smaxdum);

X = CFD_data{1,t}.x;
p_num = CFD_data{1,t}.p;
r_num = CFD_data{1,t}.rho;
T_num = CFD_data{1,t}.T;
u_num = CFD_data{1,t}.u;
M_num = CFD_data{1,t}.M;
S_num = CFD_data{1,t}.S;

X_th = theo_data{1,t}.x;
p_th = theo_data{1,t}.p;
r_th = theo_data{1,t}.rho;
T_th = theo_data{1,t}.T;
u_th = theo_data{1,t}.u;
M_th = theo_data{1,t}.M;
S_th = theo_data{1,t}.S;

% Create figures and axes
str = 'p(x) at t = '+string(time(t));
f_p = figure('Name',str,'Visible','on');
% Create toggle button !! cambia e usa radiobutton con Grid on e Grid off
sld_time_p = uicontrol('Style', 'slider',...
        'Min',timemin,'Max',timemax,'Value',time(t),...
        'Position', [400 3 120 20],...
        'Callback', @sldtime_p);     
f_p.Visible = 'on';    
%plot(X,p_num,'-o');
plot(X,p_num,'-o',X_th,p_th);
axis ([0 1 pmin-pmax*.1 pmax*(1.1)]);
xticks (0:.1:1)
xlabel('x')
ylabel('p')

str = 'rho(x) at t = '+string(time(t));
f_r = figure('Name',str,'Visible','on');
% Create toggle button !! cambia e usa radiobutton con Grid on e Grid off
sld_time_r = uicontrol('Style', 'slider',...
        'Min',timemin,'Max',timemax,'Value',time(t),...
        'Position', [400 3 120 20],...
        'Callback', @sldtime_r);     
f_r.Visible = 'on';
plot(X,r_num,'-o',X_th,r_th);
axis ([0 1 rmin-rmax*.1 rmax*(1.1)]);
xticks (0:.1:1)
xlabel('x')
ylabel('rho')

str = 'T(x) at t = '+string(time(t));
f_T = figure('Name',str,'Visible','on');
% Create toggle button !! cambia e usa radiobutton con Grid on e Grid off
sld_time_T = uicontrol('Style', 'slider',...
        'Min',timemin,'Max',timemax,'Value',time(t),...
        'Position', [400 3 120 20],...
        'Callback', @sldtime_T);     
f_T.Visible = 'on';
plot(X,T_num,'-o',X_th,T_th);
axis ([0 1 Tmin-Tmax*.1 Tmax*(1.1)]);
xticks (0:.1:1)
xlabel('x')
ylabel('T')

str = 'u(x) at t = '+string(time(t));
f_u = figure('Name',str,'Visible','on');
% Create toggle button !! cambia e usa radiobutton con Grid on e Grid off
sld_time_u = uicontrol('Style', 'slider',...
        'Min',timemin,'Max',timemax,'Value',time(t),...
        'Position', [400 3 120 20],...
        'Callback', @sldtime_u);     
f_u.Visible = 'on';
plot(X,u_num,'-o');
%plot(X,u_num,'-o',X_th,u_th);
axis ([0 1 umin-umax*.1 umax*(1.1)]);
xticks (0:.1:1)
xlabel('x')
ylabel('u')

str = 'Mach(x) at t = '+string(time(t));
f_M = figure('Name',str,'Visible','on');
% Create toggle button !! cambia e usa radiobutton con Grid on e Grid off
sld_time_M = uicontrol('Style', 'slider',...
        'Min',timemin,'Max',timemax,'Value',time(t),...
        'Position', [400 3 120 20],...
        'Callback', @sldtime_M);     
f_M.Visible = 'on';    
plot(X,M_num,'-o',X_th,M_th);
axis ([0 1 Mmin-Mmax*.1 Mmax*(1.1)]);
xticks (0:.1:1)
xlabel('x')
ylabel('M')

str = 'S(x) at t = '+string(time(t));
f_S = figure('Name',str,'Visible','on');
% Create toggle button !! cambia e usa radiobutton con Grid on e Grid off
sld_time_S = uicontrol('Style', 'slider',...
        'Min',timemin,'Max',timemax,'Value',time(t),...
        'Position', [400 3 120 20],...
        'Callback', @sldtime_S);     
f_S.Visible = 'on';    
plot(X,S_num,'-o',X_th,S_th);
axis ([0 1 Smin-Smax*.1 Smax*(1.1)]);
xticks (0:.1:1)
xlabel('x')
ylabel('S')

function sldtime_p(objtime,evntime)
   objtime.SliderStep = [dtmin/(timemax-timemin),dtmax/(timemax-timemin)];
   time_choice = objtime.Value;
   [T]=find(time<time_choice);
   t=max(T);
   T= [t,t+1];
   if (t+1>tmax)
       T=[t,t];
   end
   if (time_choice<=timemin)
       T=[1,1];
   end
   if (abs(time(T(1,1))-time_choice)>abs(time(T(1,2))-time_choice))
       t= T(2);
   else
       t= T(1);
   end
   
   timenear=time(t);
   p_num = CFD_data{1,t}.p;
   X_th = theo_data{1,t}.x;
   p_th = theo_data{1,t}.p;
   plot(X,p_num,'-o',X_th,p_th);
%   plot(X,p_num,'-o');
   str = 'p(x) at t = '+string(timenear);
   f_p.Name = str;
   axis ([0 1 pmin-pmax*.1 pmax*(1.1)]);
   xticks (0:.1:1)
   xlabel('x')
   ylabel('p')
end

function sldtime_r(objtime,evntime)
   objtime.SliderStep = [dtmin/(timemax-timemin),dtmax/(timemax-timemin)];
   time_choice = objtime.Value;
   [T]=find(time<time_choice);
   t=max(T);
   T= [t,t+1];
   if (t+1>tmax)
       T=[t,t];
   end
   if (time_choice<=timemin)
       T=[1,1];
   end
   if (abs(time(T(1,1))-time_choice)>abs(time(T(1,2))-time_choice))
       t= T(2);
   else
       t= T(1);
   end
   
   timenear=time(t);
   r_num = CFD_data{1,t}.rho;
   X_th = theo_data{1,t}.x;
   r_th = theo_data{1,t}.rho;
   plot(X,r_num,'-o',X_th,r_th);
   str = 'rho(x) at t = '+string(timenear);
   f_r.Name = str;
   axis ([0 1 rmin-rmax*.1 rmax*(1.1)]);
   xticks (0:.1:1)
   xlabel('x')
   ylabel('rho')
end

function sldtime_T(objtime,evntime)
   objtime.SliderStep = [dtmin/(timemax-timemin),dtmax/(timemax-timemin)];
   time_choice = objtime.Value;
   [T]=find(time<time_choice);
   t=max(T);
   T= [t,t+1];
   if (t+1>tmax)
       T=[t,t];
   end
   if (time_choice<=timemin)
       T=[1,1];
   end
   if (abs(time(T(1,1))-time_choice)>abs(time(T(1,2))-time_choice))
       t= T(2);
   else
       t= T(1);
   end
   
   timenear=time(t);
   T_num = CFD_data{1,t}.T;
   X_th = theo_data{1,t}.x;
   T_th = theo_data{1,t}.T;
   plot(X,T_num,'-o',X_th,T_th);
   str = 'T(x) at t = '+string(timenear);
   f_T.Name = str;
   axis ([0 1 Tmin-Tmax*.1 Tmax*(1.1)]);
   xticks (0:.1:1)
   xlabel('x')
   ylabel('T')
end

function sldtime_u(objtime,evntime)
   objtime.SliderStep = [dtmin/(timemax-timemin),dtmax/(timemax-timemin)];
   time_choice = objtime.Value;
   [T]=find(time<time_choice);
   t=max(T);
   T= [t,t+1];
   if (t+1>tmax)
       T=[t,t];
   end
   if (time_choice<=timemin)
       T=[1,1];
   end
   if (abs(time(T(1,1))-time_choice)>abs(time(T(1,2))-time_choice))
       t= T(2);
   else
       t= T(1);
   end
   
   timenear=time(t);
   u_num = CFD_data{1,t}.u;
   X_th = theo_data{1,t}.x;
   u_th = theo_data{1,t}.u;
   plot(X,u_num,'-o',X_th,u_th);
   str = 'u(x) at t = '+string(timenear);
   f_u.Name = str;
   axis ([0 1 umin-umax*.1 umax*(1.1)]);
   xticks (0:.1:1)
   xlabel('x')
   ylabel('u')
end

function sldtime_M(objtime,evntime)
   objtime.SliderStep = [dtmin/(timemax-timemin),dtmax/(timemax-timemin)];
   time_choice = objtime.Value;
   [T]=find(time<time_choice);
   t=max(T);
   T= [t,t+1];
   if (t+1>tmax)
       T=[t,t];
   end
   if (time_choice<=timemin)
       T=[1,1];
   end
   if (abs(time(T(1,1))-time_choice)>abs(time(T(1,2))-time_choice))
       t= T(2);
   else
       t= T(1);
   end
   
   timenear=time(t);
   M_num = CFD_data{1,t}.M;
   X_th = theo_data{1,t}.x;
   M_th = theo_data{1,t}.M;
   plot(X,M_num,'-o',X_th,M_th);
   str = 'Mach(x) at t = '+string(timenear);
   f_M.Name = str;
   axis ([0 1 Mmin-Mmax*.1 Mmax*(1.1)]);
   xticks (0:.1:1)
   xlabel('x')
   ylabel('M')
end

function sldtime_S(objtime,evntime)
   objtime.SliderStep = [dtmin/(timemax-timemin),dtmax/(timemax-timemin)];
   time_choice = objtime.Value;
   [T]=find(time<time_choice);
   t=max(T);
   T= [t,t+1];
   if (t+1>tmax)
       T=[t,t];
   end
   if (time_choice<=timemin)
       T=[1,1];
   end
   if (abs(time(T(1,1))-time_choice)>abs(time(T(1,2))-time_choice))
       t= T(2);
   else
       t= T(1);
   end
   
   timenear=time(t);
   S_num = CFD_data{1,t}.S;
   X_th = theo_data{1,t}.x;
   S_th = theo_data{1,t}.S;
   plot(X,S_num,'-o',X_th,S_th);
   str = 'S(x) at t = '+string(timenear);
   f_S.Name = str;
   axis ([0 1 Smin-Smax*.1 Smax*(1.1)]);
   xticks (0:.1:1)
   xlabel('x')
   ylabel('S')
end

end