close all; clear all; clc

% STM simulation
k1  = 0.031;  % k(s-2) 
k2  = 0.062;  % k(s-1) 
k3  = 0.028;  % k(s)  
k4  = 0.100;  % k(s+1)
k5  = 0.025;  % k(s+2)
kramp=0.35;   % kramp

alpha1=0.0167; alpha2=alpha1; alpha3=alpha1; alpha4=alpha1; alpha5=alpha1; 
 
p1=[3   80  120  40   4100];   %s-2
p2=[3   63  110  44   4500];   %s-1
p3=[3   60  110  54   1100];   %s
p4=[3   70  120  46   2600];   %s+1
p5=[3   60  110  52   3200];   %s+2
p6=[3  100  110  65   2600];   %boundary
p7=[2   50   70  25   3600];   %ramp

%boundary conditions
[t, q, v, o]=m27load(26, 10, 1, 1); r=occ2den(o, 4, 2, 3, 1); 
q_in=[0.75*q(:,1)]; q_in=irwsm(q_in,1,0.001); q_in=resample(q_in,6,1); 
v6=irwsm(v(:,7),1,0.001); v6=resample(v6,6,1);
v6=[110*ones(100,1); 0.95*v6(2000:4250)];
q_in=[10*ones(100,1); q_in(2160:4320)]; 
% q_urban=[2*ones(400,1); linspace(2,21,500)'; 21*ones(100,1); linspace(21,2,500)'; 2*ones(861,1);];
q_urban=[q(:,2)]; q_urban=irwsm(q_urban,1,0.001);  q_urban=1.4*resample(q_urban,6,1); q_urban=[2*ones(100,1); q_urban(2160:4320)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise other variables
r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
N_ramp=zeros(size(r1)); A=zeros(size(r3)); B=zeros(size(r3)); C=zeros(size(r3)); r3sdp_ramp=zeros(size(r3));

% -------------- Section 3 Conrol Implementation ----------------------------------------------------------
uc=[zeros(600,1); (p3(4)+5)*ones(300,1); (p3(4)-5)*ones(300,1); (p3(4)+8)*ones(300,1); p3(4)*ones(500,1)];
start=500; finish=2000;figure(3); subplot(211); plot(uc,'r'); horizon=start:6:finish; time=length(horizon);
%----------------------------------------------------------------------------------------------------------

load MCS_par %par_upstr par_jun par_down

%  for reference
%  aa=[-at(2) -at_2(2) -at_3(2)]
%  bb=[bt(1, 2) bt_2(1, 2) bt_3(1, 2)]
%  cc=[0 bt_2(2, 2) 0]
%  dd=[0 bt_2(3,2) 0]

figure(3); subplot(211); hold; subplot(212); hold; 

h2 = WAITBAR(0,'Wait MCS in progress ...')

for j=1:1000 % Number of MCS

 aa=[-par_upstr(j,1) -par_jun(j,1) -par_down(j,1)];
 bb=[par_upstr(j, 2) par_jun(j, 2) par_down(j, 2)];
 cc=[0 par_jun(j, 3) 0];
 dd=[0 par_jun(j,4) 0];
 

%Formulation of the NMSS matrix
F = [aa(1)  bb(1)      0   0;...
     bb(2)  aa(2)   cc(2)  0;...
        0   bb(3)   aa(3)  0;...
    -bb(2) -aa(2)  -cc(2)  1;];
g = [0 dd(2) 0 -dd(2)]'; 

Q=eye(size(F)); Q(4,4)=10000; Q(3,3)=1; Q(2,2)=1; Q(1,1)=1; 
r=1; k_llm = dlqri(F,g,Q,r); k_llm(end)=-k_llm(end); 
count=[];

h=waitbar(0,'Please wait.... Control in Progress (LLM) 1/7');
i=1;
for n=3:length(r2);   
    
    %link 1 [s-2]
    q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
    r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
    v1(n)=ftdcurve(p1,r1(n),1); 
    
    %link 2 [s-1]
    q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
    r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
    v2(n)=ftdcurve(p2,r2(n),1);
    
    %Ramp   
    q_ramp(n)=alpha1*min(v_ramp(n-1),v3(n-1))*r_ramp(n-1);
    r_ramp(n)=r_ramp(n-1)+((q_urban(n-1)-q_ramp(n-1))*kramp);    
    v_ramp(n)=ftdcurve(p7,r_ramp(n),1);
    N_ramp(n)=N_ramp(n-1)+alpha1*(q_urban(n-1)-q_ramp(n-1));
       
    %link 3 [s]
    q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
    r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
    v3(n)=ftdcurve(p3,r3(n),1);
    
    %link 4 [s+1]
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);if (n>1620)&(n<1680); q4(n)=1/3*q4(n); end
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    v4(n)=ftdcurve(p4,r4(n),1);
    
    %link 5 [s+2]
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
        
    if  (n>=horizon(1)) & (n<horizon(end))      
        if (i<time+1)&(n==horizon(i))    
            i=i+1;
            
            e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + k_llm(4)*e(n) - k_llm(1)*(r2(n)-r2(n-1)) - k_llm(2)*(r3(n)-r3(n-1)) - k_llm(3)*(r4(n)-r4(n-1));
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
        q_ramp(n)=q_urban(n);
  
        waitbar(n/length(q1),h);
    end
end   
close(h)

if r3(1800)>uc(1800)+20; count=[count j];end
figure(3); 
subplot(211); h1=plot(r3,'k'); 
subplot(212); h1=plot(q_ramp,'k'); 

waitbar(j/1000,h2)
end
close(h2)

subplot(211);set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['05:00'; '06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'])  
h3=ylabel('Density (veh/km)'); set(h3,'FontSize',14);h4=xlabel('Time (hr)'); set(h4,'FontSize',14);grid
axis([0 2200 0 80]);h=plot(uc,'r'); set(h,'LineWidth',2); h=plot(p3(4)*ones(2100,1),'-.'); set(h,'LineWidth',2);
subplot(212);set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['05:00'; '06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'])  
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);axis([0 1200 0 30]); h3=ylabel('Flow (veh/min)'); set(h3,'FontSize',14);
axis([0 2200 0 30]);grid
% saveas(gcf, 'output', 'bmp')
