%% Srijan Dasgupta - Assignment 2a
clc
close all
clear all

%% 1D ground problem 

%% 1.Input Data 
% Physical Data 
L=2;                                            %The depth of the ground(m) 
x1=0;                                           %length start           (m)
x2=L;                                           %length end             (m)
H=1;                                            %unit height            (m)
W=1;                                            %unit width             (m) 
ro=2300;                                        %density                (kg/m3)
cp=700;                                         %specific heat          (J/kgK)
lamda=1.6;                                      %thermal conductivity   (W/mK)
alpha=8.5;                                      %heat coefficient       (W/m2K)
K=273.15;                                       %degree to kelvin conversion 
T0=19+K;                                        %initial temperature at t=0  (K)
Tw=19+K;                                        %isothermal wall temperature (K)
A0=19;                                          %Coefficient    (degree celcius)
A1=5.7;                                         %Coefficient    (degree celcius)
A2=9.1;                                         %Coefficient    (degree celcius)
w1=(2*pi)/(24*3600);                            
w2=(2*pi)/(365*24*3600);
Text= @(x) A0+A1*sin(w1*x)+A2*sin(w2*x);        %external air temperature function 
% Define time 
tstart=0;                                       %start time
tend=100*24*3600;                               %end time 
delt=10000;                                     %timestep 
% Numerical Data 
N=101;                                          %number of CVs
beta=1.0;                                       %transient solution coefficient
tolerance= 1e-9;                                %tolerance of error 
maxiter=1e6;                                    %iteration limit 
%% 2.Previous Calculations 
% Domain discretization 
x_cv=linspace(x1,x2,N+1);                       %Face position matrix    (m)
x_p=zeros(1,N+2);                               %node postion matrix     (m)
V=zeros(1,N+2);                                 %volume matrix for CVs   (m3) 
S=ones(1,N+1);                                  %cross sectional area    (m2)
x_p(1)=x1;                                      %1st node position       (m)
x_p(end)=x2;                                    %last node position      (m)
del_x=(x2-x1)/N;                                %discretization seperation thickness (m)

for i=2:1:N+1
    x_p(i)=(x_cv(i)+x_cv(i-1))/2;               %internal node positions (m)
    V(i)=(x_cv(i)-x_cv(i-1))*W*H;               %Volume of CVs
end 

% Define coefficient vectors 
ap=ones(1,N+2);
ae=zeros(1,N+2);
aw=zeros(1,N+2);
bp=zeros(1,N+2);

%% 3.Initial Temperature 
T0=T0*ones(1,N+2);                              %initial temperature (K)
Tn=T0;                                          %current timestep temperature 
%% 4.Calculation Phase  
tn=tstart;                                      %current time 
tplus=tn+delt;                                  %t(n+1)=tn+delt=next timestep

% 4.1 Estimated temperature at current instant
Tplus=Tn;                                       %next timestep guess temperature 
Tguess=Tplus;                                   %store Tplus for iterative comparison

%Time loop 
j=1; 
while tplus<=tend 
    % 4.2 Estimation of the coefficient values using Tn and Tn+1
    time(j,1)=tplus;                                    %storing time
    Tplus(j,:)=Tn;                                      %guess temperature of next timestep 
    %coefficients for first nodes 
    aw(1)=0;
    ae(1)=lamda/(x_p(2)-x_p(1));
    ap(1)=ae(1)+alpha; 
    bp(1)=alpha*(Text(tplus)+K);
    %coefficients for internal nodes 
    for i=2:1:N+1
        aw(i)=(lamda*S(i-1))/(x_p(i)-x_p(i-1));
        ae(i)=(lamda*S(i))/(x_p(i+1)-x_p(i));
        ap(i)=beta*ae(i)+beta*aw(i)+((ro*V(i)...
            *cp)/delt);
        bp(i)=((ro*V(i)*cp)/delt)*Tn(i)+(1-beta)...
            *(-(ae(i)+aw(i))*Tn(i)+aw(i)*Tn(i-1)...
            +ae(i)*Tn(i+1));
    end
    ap(end)=1;
    aw(end)=0;
    ae(end)=0;
    bp(end)=Tw;
    % Node temperature loop
    iter=0; 
    residual=1;
    % 4.3 Solve the discretization equations 
    while iter<maxiter                                          
       Tplus(j,1)=(ae(1)*Tplus(j,2)+bp(1))/ap(1);
       Tplus(j,end)=bp(end);
       for i=2:1:N+1
           Tplus(j,i)=(beta*ae(i)*Tplus(j,i+1)+beta*aw(i)...
               *Tplus(j,i-1)+bp(i))/ap(i);
       end
       
       % 4.4 Is max|Tn+1-Tn+1*|<tolerance 
       residual=max(abs(Tplus(j,:)-Tguess));        %checking residual after iteration
       if residual>tolerance 
           Tguess=Tplus(j,:);                       % Feedback=No, Tn+1*=Tn+1
       else 
           break;                                   %Feedback=Yes, Go to next time step 
       end 
       iter=iter+1; 
    end
    figure(1)
    plot(x_p,Tplus(j,:));                           %plotting temperature for each timestep
    ylabel('Temperature') % left y-axis 
    xlabel('distance(m)')
    title('Temperature Variation (transient case)');
    drawnow
    %% 5. New time instant?
    tplus=tplus+delt;                               % Feedback=yes
    Tn=Tplus(j,:);                                  % Feedback=yes
    j=j+1;
end 
%% Energy Balance (for last two timesteps)
%Heat accumulation inside all the nodes between time tn+1 and tn
for k=2:1:length(x_p)-1                             %internal nodes 
    Qacc(k)=ro*cp*V(k)*((Tplus(end,k)-Tplus...      
        (end-1,k))/(time(end)-time(end-1)));
end
Qacc=sum(Qacc);                                     %accumulation 

%Heat transfer from air to the first node for last timestep 
if beta==1
    Qe=alpha*(Text(time(end))+K-Tplus(end,1))*S(1);
elseif beta==0.5
    Te_avg=(Text(time(end))+K+Text(time(end-1))+K)/2;
    Tplus_avg=(Tplus(end,1)+Tplus(end-1,1))/2;
    Qe=alpha*(Te_avg-Tplus_avg)*S(1);
end 

%Heat coming from isothermal point at x=L for last timestep
if beta==1
    Qc=-lamda*((Tplus(end,end-1)-Tplus(end,end))/(x_p(end)-x_p(end-1)));
elseif beta==0.5
    Tplus1_avg=(Tplus(end,end-1)+Tplus(end-1,end-1))/2;
    Qc=-lamda*((Tplus1_avg-Tplus(end,end))/(x_p(end)-x_p(end-1)));
end 

% Accumulation = Heat from air + heat from isothermal wall
Energy_balance=Qacc-Qe-Qc;                          %Local time energy balance 
%% Variation of Air temperature with time 
j=1;
k=0:delt:100*24*3600;                               %time evolution 
for i=0:delt:100*24*3600
     Te(j)=Text(i);                                 %external air temperature
     j=j+1;
end 
figure (2)                                          %plotting air temperature vs time 
plot (k,Te)
xlabel('time (sec)')
ylabel('Temperature (Celcius)')
title ('Evolution of Air temperature over time')
%% Energy Balance for t=0 to t=100*24*3600 (Implicit Method)
n=((tend-tstart)/delt)+1;
Te_avg=((sum(Te))/(n))+K;                           %average temeperature of the air 
Tplus_1avg=(sum(Tplus(:,1))+T0(1))/(n);             %average of the first node temperature
Tp_end1_avg=(sum(Tplus(:,end-1))+T0(end-1))/(n);    %average temperature of the 2nd last node from right

%heat from air 
Qe_total=alpha*(Te_avg-Tplus_1avg);                 %external air boundary
%heat from isothermal boundary
Qc_total=-lamda*(Tp_end1_avg-Tplus(end,end))...     %isothermal boundary
    /(x_p(end)-x_p(end-1));

%accumulation 
for k=2:1:length(x_p)-1                             
    Qacc_total(k)=ro*cp*V(k)*((Tplus(end,k)-T0(k)));
end
Qacc_total=sum(Qacc_total)/(tend-tstart);           %total accumulation 
EB_total=Qacc_total-Qc_total-Qe_total;              %global time energy balance 
%% Plots 
figure (3)                                          %final temperature distribution 
plot (x_p,Tplus(end,:))
xlabel('Distance (m)')
ylabel('Temperature (K)')
title ('Nodal temperatures at final timestep')

figure (4)                                           %Temperature at x=0
plot (time,Tplus(:,1))
xlabel('time (sec)')
ylabel('Temperature (K)')
title ('Temperature evolution at x=0')

%calculating temperature at x=L/2
if rem(N,2)==0
    n1=(N+2)/2;
    n2=(N+4)/2; 
    Tn1=Tplus(:,n1);
    Tn2=Tplus(:,n2);
    TL_2=(Tn1+Tn2)/2;
else 
    n=(N+3)/2;
    TL_2=Tplus(:,n);     
end                     
figure (5)                                          %Temperature at x=L/2
plot (time,TL_2)
xlabel('Time (sec)')
ylabel('Temperature (K)')
title ('Temperature evolution at x=L/2')


figure (6)                                          %isothermal node temperature evolution 
plot (time,Tplus(:,end))
xlabel('Time (sec)')
ylabel('Temperature (K)')
title ('Temperature evolution at x=L')

