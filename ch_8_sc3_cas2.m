close all; clear all; clc


% STM simulation
k1  = 0.031;  % k(s-2) 
k2  = 0.062;  % k(s-1) 
k3  = 0.028;  % k(s)  
k4  = 0.100;  % k(s+1)
k5  = 0.025;  % k(s+2)
kramp=0.35;   % kramp

alpha1=0.0167; alpha2=alpha1; alpha3=alpha1; alpha4=alpha1; alpha5=alpha1; 

p1=[3   90  120  40   1100];   %s-2
p2=[3   63  110  42   3000];   %s-1
p3=[3   65  110  45   3100];   %s
p4=[3   70  110  45   1600];   %s+1
p5=[3   80  110  35   1600];   %s+2
p7=[2   40   60  25   3600];   %ramp


%boundary conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t, q, v, o]=m27load(13, 11, 1, 1); r=occ2den(o, 4, 2, 3, 1); 
q_in=[0.55*q(:,1)]; q_in=irwsm(q_in,1,0.001); q_in=resample(q_in,6,1); 
v6=irwsm(v(:,7),1,0.001); v6=resample(v6,6,1);
v6=[110*ones(100,1); v6(2000:4000)];
q_in=[10*ones(100,1); q_in(2000:4000)]; 
q_urban=[q(:,2)]; q_urban=irwsm(q_urban,1,0.001);  q_urban=1.2*resample(q_urban,6,1); 
q_urban=[2*ones(100,1); q_urban(2000:4000)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initialise other variables
r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
N_ramp=zeros(size(r1));
A=zeros(size(r3)); B=zeros(size(r3)); C=zeros(size(r3));
r3sdp_ramp=zeros(size(r3));


% === Open Loop Simulation ===
h=waitbar(0,'Open Loop Simulation ');

for n=3:length(q_in);
    
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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1); 
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4); 
    v4(n)=ftdcurve(p4,r4(n),1);     
    
    %link 5 [s+2]
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
    
    waitbar(n/length(r2),h);
end

close(h)

figure(3);subplot(211);hold;grid
h1=plot(r3,'k');                         set(h1,'LineWidth',2)
h2=plot(p3(4)*ones(length(q1),1),'k:');  set(h2,'LineWidth',2)
h3=ylabel('Density (veh/km)');           set(h3,'FontSize',14);
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
set(gca, 'xtick',[1 380 760 1140 1300 1520 1900], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '09:27'; '10:00'; '11:00';])  

subplot(212); hold;grid
plot(q_urban,'k'); plot(q_ramp,'k');
set(gca, 'xtick',[1 380 760 1140 1300 1520 1900], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '09:27'; '10:00'; '11:00';])  
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
h3=ylabel('Flow (veh/min)'); set(h3,'FontSize',14);
figure(6);
subplot(321); plot(r1,q1,'.'); hold; plot(p1(4),linspace(0,120),'k.'); plot(p1(4),max(q1),'r.');title('s-2');
subplot(322); plot(r2,q2,'.'); hold; plot(p2(4),linspace(0,120),'k.'); plot(p2(4),max(q2),'r.');title('s-1 {load}');
subplot(323); plot(r3,q3,'.'); hold; plot(p3(4),linspace(0,120),'k.'); plot(p3(4),max(q3),'r.');title('s {junction}');
subplot(324); plot(r4,q4,'.'); hold; plot(p4(4),linspace(0,120),'k.'); plot(p4(4),max(q4),'r.');title('s+1');
subplot(325); plot(r5,q5,'.'); hold; plot(p5(4),linspace(0,120),'k.'); plot(p5(4),max(q5),'r.');title('s+2');
subplot(326); plot(r_ramp,q_ramp,'.'); hold; plot(p7(4),linspace(0,120),'k.'); plot(p7(4),max(q_ramp),'r.');title('ramp');
xlabel('Density veh/km'); ylabel('Flow veh/h');

%time delays
r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r1)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4)); r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5));
v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v4)); v4_2d = lag(v4,2,mean(v4)); v5_2d = lag(v5,2,mean(v5));


%Period density is above critical
TTT_open = (sum(q_ramp(835:1402)+q2(835:1402)-q5(835:1402)) + q2(834)+q_ramp(834)-q5(834) ) * (1402-834)/3600;
TWT_open = sum((q_urban(780:1128)-q_ramp(780:1128))) * (1128-780)/3600;
TTS_open=TTT_open+TWT_open;

% -------------- Section 3 Conrol Implementation -------------
uc=[zeros(200,1); (p3(4))*ones(2000,1);];  figure(3); subplot(211); plot(uc,'r'); 
horizon=876:1300; time=length(horizon);
%------------------------------------------------------------------

% *** LLM Controller on STM model ***
load llm_par

%Formulation of the NMSS matrix
F = [aa(1)  bb(1)      0   0;...
     bb(2)  aa(2)   cc(2)  0;...
        0   bb(3)   aa(3)  0;...
    -bb(2) -aa(2)  -cc(2)  1;];
g = [0 dd(2) 0 -dd(2)]'; 

Q=eye(size(F)); Q(4,4)=1000; Q(3,3)=1; Q(2,2)=1; Q(1,1)=1; 
r=1; k_llm = dlqri(F,g,Q,r); k_llm(end)=-k_llm(end); 

% Qb=1; C=[1 1 100 1000000]; 
% Q=C'*Qb*C; 
% r=1; k_llm = dlqri(F,g,Q,r); k_llm(end)=-k_llm(end); 
% 


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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
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
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7;  end  
            
%             if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
%         if N_ramp(n)>93; q_ramp(n)=14.7; 
%         elseif (N_ramp(n)<=93); q_ramp(n)=0.5*N_ramp(n); 
%             if q_ramp(n)>16.7; q_ramp(n)=14.7; end   
%         elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=14.7; end 
%         end
q_ramp(n)=min((40/r3(n))*16.7,16.7);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end   

    end
    waitbar(n/length(q1),h);    
end
close(h)

figure(3); 
subplot(211); h1=plot(r3,'k'); h5=legend('Open loop Density','Critical Density','Reference Level','Closed loop Density');
plot(1300,linspace(0,150,'r'))
subplot(212); h1=plot(q_ramp,'k'); h5=legend('Urban Demand', 'Unmetered Ramp Flow', 'Metered Ramp Flow');

% TRAFFIC Criteria
TTT_llm = (sum(q_ramp(876:1300)+q2(876:1300)-q5(876:1300)) + q2(875)+q_ramp(875)-q5(875) ) * (1300-876)/3600;
TWT_llm = sum((q_urban(876:1300)-q_ramp(876:1300))) * (1300-876)/3600;
TTS_llm=TTT_llm+TWT_llm;

TTT_llm_add = (sum(q_ramp(1301:1608)+q2(1301:1608)-q5(1301:1608)) + q2(1300)+q_ramp(1300)-q5(1300) ) * (1608-1301)/3600;



% % Control Criteria
IAE_llm=sum(abs(r3(876:1300)-uc(876:1300)));
resid=r3(876:1300)-uc(876:1300); index=find(resid>0); IAE2_llm=sum(resid(index));
load_llm=sum(q_ramp(876:1300))/6
index=find(r3(876:1300)<=uc(end)+0.05);Set_time_llm=index(1);


%second order llm control
load llm2_par

% Formulation of the NMSS matrices
F_llm2=[ 1   aa2(1)  1  cc2(1)  0     0       0 ; ...
        1    0      0   0      0     0       0 ; ...     
        0   bb2(2)  1  aa2(2)  1   cc2(2)    0 ; ...
        0     0     1   0      0     0       0 ; ...
        0     0     0  bb2(3)  1   aa2(3)    0 ; ...
        0     0     0   0      1     0       0 ; ...
        0  -bb2(2) -1 -aa2(2) -1  -cc2(2)    1 ];

g_llm2 =[0 0 dd2(1) 0 0 0 -dd2(1)]';

Q=eye(size(F_llm2)); Q(end,end)=1000;
Q(1,1)=1; Q(2,2)=1; Q(3,3)=1; Q(4,4)=1; Q(5,5)=1; Q(6,6)=1;

r=1; k_llm2 = dlqri(F_llm2,g_llm2,Q,r); k_llm2(end)=-k_llm2(end); 

% Qb=1; C=[1 1 1 1 1000 100 10000]; 
% Q=C'*Qb*C; Q=Q+eye(7);
% r=1; k_llm2 = dlqri(F_llm2,g_llm2,Q,r); k_llm2(end)=-k_llm2(end); 



% *** LLM Controller on STM model ***

h=waitbar(0,'Please wait.... Control in Progress (LLM 7 states) 2/7');
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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    v4(n)=ftdcurve(p4,r4(n),1);
    
    %link 5 [s+2]
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
    
    if  (n>=horizon(1)) & (n<horizon(end))      
        if (i<time+1)&(n==horizon(i))    
            i=i+1;
            %81 times
            e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + k_llm2(7)*e(n) - k_llm2(1)*(r2(n)-r2(n-1)) ...
                - k_llm2(2)*(r2(n-1)-r2(n-2)) ...
                - k_llm2(3)*(r3(n)-r3(n-1)) ...
                - k_llm2(4)*(r3(n-1)-r3(n-2)) ...
                - k_llm2(5)*(r4(n)-r4(n-1)) ...
                - k_llm2(6)*(r4(n-1)-r4(n-2));  
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
            if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
%         if N_ramp(n)>93; q_ramp(n)=14.7; 
%         elseif (N_ramp(n)<=93); q_ramp(n)=0.5*N_ramp(n); 
%             if q_ramp(n)>16.7; q_ramp(n)=14.7; end   
%         elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=14.7; end 
%         end
q_ramp(n)=min((40/r3(n))*16.7,16.7);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end   

    end
    waitbar(n/length(q1),h);    
end
close(h)

figure(3); 
subplot(211);h1=plot(r3); 
subplot(212); plot(q_ramp); 
% TRAFFIC Criteria
TTT_llm2 = (sum(q_ramp(876:1300)+q2(876:1300)-q5(876:1300)) + q2(875)+q_ramp(875)-q5(875) ) * (1300-876)/3600;
TWT_llm2 = sum((q_urban(876:1300)-q_ramp(876:1300))) * (1300-876)/3600;
TTS_llm2=TTT_llm2+TWT_llm2;

% % Control Criteria
IAE_llm2=sum(abs(r3(876:1300)-uc(876:1300)));
resid=r3(876:1300)-uc(876:1300); index=find(resid>0); IAE2_llm2=sum(resid(index));
load_llm2=sum(abs(q_ramp(876:1300))/6);
index=find(r3(876:1300)<=uc(end)+0.05);Set_time_llm2=index(1);

% *** LLM Controller on STM model ***

h=waitbar(0,'Please wait.... Control in Progress (alinea 4/7)');
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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
    v4(n)=ftdcurve(p4,r4(n),1);
    
    %link 5 [s+2]
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
    
    if  (n>=horizon(1)) & (n<horizon(end))      
        if (i<time+1)&(n==horizon(i))    
            i=i+1;
            %81 times
            e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + 13*e(n);
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
            if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
%         if N_ramp(n)>93; q_ramp(n)=16.7; 
%         elseif (N_ramp(n)<=93); q_ramp(n)=0.5*N_ramp(n); 
%             if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
%         elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
%         end
q_ramp(n)=min((40/r3(n))*16.7,16.7);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end   

    end
    waitbar(n/length(q1),h);
    
end
close(h)

% figure(3); subplot(211);h1=plot(r3,'m'); legend('Open Loop Density', 'Critical Density','Reference Level','LLM-PIP','ALINEA');
% subplot(212);plot(q_ramp,'m'); legend('Urban Demand','Unmetered Ramp Output','LLM-PIP','ALINEA');


% TRAFFIC Criteria
TTT_alinea = (sum(q_ramp(876:1300)+q2(876:1300)-q5(876:1300)) + q2(875)+q_ramp(875)-q5(875) ) * (1300-876)/3600;
TWT_alinea = sum((q_urban(876:1300)-q_ramp(876:1300))) * (1300-876)/3600;
TTS_alinea=TTT_alinea+TWT_alinea;

TTT_alinea_Add = (sum(q_ramp(1300:1800)+q2(1300:1800)-q5(1300:1800)) + q2(1299)+q_ramp(1299)-q5(1299) ) * (1800-1300)/3600;


% % Control Criteria
IAE_alinea=sum(abs(r3(876:1300)-uc(876:1300)));
resid=r3(876:1300)-uc(876:1300); index=find(resid>0); IAE2_alinea=sum(resid(index));
load_alinea=sum(abs(q_ramp(876:1300))/6);
index=find(r3(876:1300)<=uc(end)+0.05);Set_time_alinea=index(1);


ALINEA=['    [1]  ALINEA       '  num2str(IAE_alinea)  '     ' num2str(IAE2_alinea)  '     '  num2str(load_alinea) ...
        '       ' num2str(load_alinea) '               ' num2str(Set_time_alinea)];
LLM   =['    [2]  LLM          '  num2str(IAE_llm)  '      ' num2str(IAE2_llm)  '    '  num2str(load_llm) ...
        '       ' num2str(load_llm) '               ' num2str(Set_time_llm)];
LLM_d =['    [3]  LLM_det      '  num2str(IAE_llm2)  '      ' num2str(IAE2_llm2)  '    '  num2str(load_llm2) ...
        '       ' num2str(load_llm2) '                ' num2str(Set_time_llm2)];

fprintf('                   Control Criteria                    Results                      \n \n');
fprintf('      Controller       IAE        Sum(r +)   Ramp Effort     Load Reduct       Settling time    \n  ');
disp('---------------------------------------------------------------------------------------------------');
disp(ALINEA)
disp(LLM)
disp(LLM_d)



ALINEA=['    [1]  ALINEA       '  num2str(TTT_alinea)  '      ' num2str(TWT_alinea)  '    '  num2str(TTS_alinea)];
LLM   =['    [2]  LLM          '  num2str(TTT_llm)  '      ' num2str(TWT_llm)  '    '  num2str(TTS_llm) ];
LLM_d =['    [3]  LLM_det      '  num2str(TTT_llm2)  '      ' num2str(TWT_llm2)  '    '  num2str(TTS_llm2)];

fprintf('                  Traffic Criteria:                    Results                      \n \n');
fprintf('      Controller       TTT           TWT          TTS    \n  ');
disp('---------------------------------------------------------------------------------------------------');
disp(ALINEA)
disp(LLM)
disp(LLM_d)
