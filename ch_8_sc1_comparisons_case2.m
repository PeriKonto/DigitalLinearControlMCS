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
p2=[3   63  110  42.5   3000];   %s-1
p3=[3   65  110  50   3100];   %s
p4=[3   90  120  35   3600];   %s+1
p5=[3   80  110  70   1600];   %s+2
p6=[3  100  110  65   2600];   %boundary
p7=[2   40   60  25   3600];   %ramp

%Boundary Conditions
%Downstream
r6=[10*ones(150,1); linspace(10,70,150)'; 70*ones(700,1); linspace(70,10,150)'; 10*ones(260,1)];  v6=ftdcurve(p6,r6,1);

%Upstream
q_in=[10*ones(150,1); linspace(10,70,150)'; 70*ones(500,1); linspace(70,10,150)'; 10*ones(260,1)]; q_in=0.70*q_in;  

%Ramp Input
q_urban=[2*ones(300,1); linspace(2,20,150)'; 20*ones(200,1); linspace(20,2,150)'; 2*ones(660,1);];

% initialise other variables
r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
N_ramp=zeros(size(r1));


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
set(gca, 'xtick',[1 380 760 1140 ], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; ])  

subplot(212); hold;grid
plot(q_urban,'k'); plot(q_ramp,'k');set(gca, 'xtick',[1 380 760 1140 ], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; ])  
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
h3=ylabel('Flow (veh/10*s)'); set(h3,'FontSize',14);


%time delays
r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r1)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4)); r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5));
v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v4)); v4_2d = lag(v4,2,mean(v4)); v5_2d = lag(v5,2,mean(v5));

%Period density is above critical
TTT_open = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(400)+q_ramp(400)-q3(400) ) * (960-400)/3600;
TWT_open = sum((q_urban(400:676)-q_ramp(400:676))) * (676-400)/3600;
TTS_open=TTT_open+TWT_open;

r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
N_ramp=zeros(size(r1)); r3_llm=NaN*ones(size(q_in)); r3_llm(1:2)=1;



% -------------- Section 3 Conrol Implementation -------------
% uc=[zeros(200,1); (p3(4)-10)*ones(400,1); (p3(4)+10)*ones(100,1); (p3(4))*ones(800,1)];  
uc=[zeros(200,1); (p3(4)-5)*ones(1000,1)];  figure(3); subplot(211); plot(uc,'r'); 
horizon=400:6:960; time=length(horizon);
%------------------------------------------------------------------


% *** LLM Controller on STM model ***
load llm_par

%Formulation of the NMSS matrix
F = [aa(1)  bb(1)      0   0;...
     bb(2)  aa(2)   cc(2)  0;...
        0   bb(3)   aa(3)  0;...
    -bb(2) -aa(2)  -cc(2)  1;];

g = [0 dd(2) 0 -dd(2)]'; 

% F = [aa(1)  bb(1)      0   0;...
%      bb(2)  aa(2)   cc(2)  0;...
%         0   bb(3)   aa(3)  0;...
%         0  -bb(3)  -aa(3)  1;];
% 
% g = [0 dd(2) 0 0]'; 

Q=eye(size(F)); Q(4,4)=1000; Q(3,3)=1; Q(2,2)=1; Q(1,1)=1; 
%Q(4,4) was 100000 but I changed to 1000 sicne this is how much I saw on the thesis
r=1; k_llm = dlqri(F,g,Q,r); k_llm(end)=-k_llm(end); 


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
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1); %if (n>500)&(n<600); r5(n-1)=1.2*r5(n-1); end
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
    
    %     r3_llm(n)=aa(2)*r3(n-1)+bb(2)*r2(n-1)+cc(2)*r4(n-1)+dd(2)*q_ramp(n-1);
    
    if  (n>=horizon(1)) & (n<horizon(end))      
        if (i<time+1)&(n==horizon(i))    
            i=i+1
           e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + k_llm(4)*e(n) - k_llm(1)*(r2(n)-r2(n-1)) - k_llm(2)*(r3(n)-r3(n-1)) - k_llm(3)*(r4(n)-r4(n-1));
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7;  end  
            
            if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
        if N_ramp(n)>40; q_ramp(n)=16.7; %Release cars at capacity level
        elseif (N_ramp(n)<=40)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
        elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
        end
    end
    waitbar(n/length(q1),h);
    
end
close(h)

figure(3); 
subplot(211); h1=plot(r3,'k'); h5=legend('Open loop Density','Critical Density','Reference Level','Closed loop Density');
subplot(212); h1=plot(q_ramp,'k'); h5=legend('Urban Demand', 'Unmetered Ramp Flow', 'Metered Ramp Flow');

% TRAFFIC Criteria
TTT_llm = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(399)+q_ramp(399)-q4(399) ) * (960-400)/3600;
TWT_llm = sum((q_urban(400:739)-q_ramp(400:739))) * (739-400)/3600;
TTS_llm=TTT_llm+TWT_llm;

% % Control Criteria
IAE_llm=sum(abs(r3(400:960)-uc(400:960)));
resid=r3(400:960)-uc(400:960); index=find(resid>0); IAE2_llm=sum(resid(index));
load_llm=sum(abs(q_ramp(400:960)));
index=find(r3(400:960)<=uc(end)+0.05);Set_time_llm=index(1);


% break
% 
% 
% 
% %second order llm control
% 
% load llm2_par
% 
% % Formulation of the NMSS matrices
% F_llm2=[ 1   aa2(1)  1  cc2(1)  0     0       0 ; ...
%         1    0      0   0      0     0       0 ; ...     
%         0   bb2(2)  1  aa2(2)  1   cc2(2)    0 ; ...
%         0     0     1   0      0     0       0 ; ...
%         0     0     0  bb2(3)  1   aa2(3)    0 ; ...
%         0     0     0   0      1     0       0 ; ...
%         0  -bb2(2) -1 -aa2(2) -1  -cc2(2)    1 ]
% 
% g_llm2 =[0 0 dd2(1) 0 0 0 -dd2(1)];
% 
% Q=eye(size(F_llm2)); Q(end,end)=1000000;
% Q(1,1)=1; Q(2,2)=1; Q(3,3)=100000; Q(4,4)=10000; Q(5,5)=1; Q(6,6)=1;
% 
% r=1; k_llm2 = dlqri(F_llm2,g_llm2',Q,r); k_llm2(end)=-k_llm2(end); 
% 
% Qb=1; C=[1 1 1000 1000 1 1 1000000]; 
% Q=C'*Qb*C; Q=Q+eye(7);
% r=1; k_llm2 = dlqri(F_llm2,g_llm2',Q,r); k_llm2(end)=-k_llm2(end); 
% 
% 
% 
% % *** LLM Controller on STM model ***
% 
% h=waitbar(0,'Please wait.... Control in Progress (LLM 7 states) 2/7');
% i=1;
% for n=3:length(r2);   
%     
%     %link 1 [s-2]
%     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
%     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
%     v1(n)=ftdcurve(p1,r1(n),1); 
%     
%     %link 2 [s-1]
%     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
%     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
%     v2(n)=ftdcurve(p2,r2(n),1);
%     
%     %Ramp   
%     q_ramp(n)=alpha1*min(v_ramp(n-1),v3(n-1))*r_ramp(n-1);
%     r_ramp(n)=r_ramp(n-1)+((q_urban(n-1)-q_ramp(n-1))*kramp);    
%     v_ramp(n)=ftdcurve(p7,r_ramp(n),1);
%     N_ramp(n)=N_ramp(n-1)+alpha1*(q_urban(n-1)-q_ramp(n-1));
%     
%     
%     %link 3 [s]
%     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
%     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
%     v3(n)=ftdcurve(p3,r3(n),1);
%     
%     %link 4 [s+1]
%     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
%     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
%     v4(n)=ftdcurve(p4,r4(n),1);
%     
%     %link 5 [s+2]
%     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
%     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
%     v5(n)=ftdcurve(p5,r5(n),1);
%     
%     if  (n>=horizon(1)) & (n<horizon(end))      
%         if (i<time+1)&(n==horizon(i))    
%             i=i+1
%             %81 times
%             e(n) = uc(n) - r3(n);
%             
%             q_ramp(n) = q_ramp(n-1) + k_llm2(7)*e(n) 
%                 - k_llm2(1)*(r3(n)-r3(n-1)) ...
%                 - k_llm2(2)*(r3(n-1)-r3(n-2)) ...
%                 - k_llm2(3)*(r2(n)-r2(n-1)) ...
%                 - k_llm2(4)*(r2(n-1)-r2(n-2)) ...
%                 - k_llm2(5)*(r4(n)-r4(n-1)) ...
%                 - k_llm2(6)*(r4(n-1)-r4(n-2));  
%             
%             if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
%             if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
%             
%         else   
%             q_ramp(n)=q_ramp(n-1);
%         end
%     end  
%     if  n>=horizon(time)  
%         if N_ramp(n)>40; q_ramp(n)=16.7; %Release cars at capacity level
%         elseif (N_ramp(n)<=40)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
%         elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
%         end
%     end
%     waitbar(n/length(q1),h);
%     
% end
% close(h)
% 
% figure(3); 
% subplot(211);h1=plot(r3); 
% subplot(212); plot(q_ramp); 
% 
% % TRAFFIC Criteria
% TTT_llm2 = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(399)+q_ramp(399)-q4(399) ) * (960-400)/3600;
% TWT_llm2 = sum((q_urban(400:739)-q_ramp(400:739))) * (739-400)/3600;
% TTS_llm2=TTT_llm2+TWT_llm2;
% 
% % % Control Criteria
% IAE_llm2=sum(abs(r3(400:960)-uc(400:960)));
% resid=r3(400:960)-uc(400:960); index=find(resid>0); IAE2_llm2=sum(resid(index));
% load_llm2=sum(abs(q_ramp(400:960)));
% index=find(r3(400:960)<=uc(end)+0.05);Set_time_llm2=index(1);
% 
% 
% % adaptive control
% r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
% v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
% q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
% r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
% v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
% v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
% N_ramp=zeros(size(r1));
% 
% P1_up=zeros(94,1); P1_us_l=zeros(94,1); P2_up=zeros(94,1); P2_us_l=zeros(94,1); P3_us_l=zeros(94,1); P3_up=zeros(94,1);
% %============== link s-1
% Q1_up=[]; P1_up(1)=10000; P1_us_l(1)=10000; Q2_up=[]; P2_up(1)=10000; P2_us_l(1)=10000; Q3_up=[]; P3_up(1)=10000; P3_us_l(1)=1000; 
% P1_j_l=zeros(94,1); P2_j_l=zeros(94,1); P3_j_l=zeros(94,1); P1_j=zeros(94,1); P2_j=zeros(94,1); P3_j=zeros(94,1); P4_j=zeros(94,1); P4_j_l=zeros(94,1);
% %============== link s
% Q1_j=[]; P1_j_l(1)=10000; P1_j(1)=10000; Q2_j=[]; P2_j_l(1)=10000;  P2_j(1)=10000;
% Q3_j=[]; P3_j_l(1)=10000; P3_j(1)=10000; Q4_j=[]; P4_j_l(1)=10000;  P4_j(1)=10000;
% P1_dw_l=zeros(94,1); P1_dw=zeros(94,1); P2_dw_l=zeros(94,1); P2_dw=zeros(94,1);
% %=============== link s+1
% Q1_dw=[]; P1_dw_l(1)=10000; P1_dw(1)=10000; Q2_dw=[]; P2_dw_l(1)=10000; P2_dw(1)=10000;
% 
% h=waitbar(0,'Please wait.... Control in adaptive (LLM) 3/7');
% i=1;
% for n=3:length(r2);   
%     
%     %link 1 [s-2]
%     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
%     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
%     v1(n)=ftdcurve(p1,r1(n),1); 
%     
%     %link 2 [s-1]
%     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
%     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
%     v2(n)=ftdcurve(p2,r2(n),1);
%     
%     %Ramp   
%     q_ramp(n)=alpha1*min(v_ramp(n-1),v3(n-1))*r_ramp(n-1);
%     r_ramp(n)=r_ramp(n-1)+((q_urban(n-1)-q_ramp(n-1))*kramp);    
%     v_ramp(n)=ftdcurve(p7,r_ramp(n),1);
%     N_ramp(n)=N_ramp(n-1)+alpha1*(q_urban(n-1)-q_ramp(n-1));
%     
%     
%     %link 3 [s]
%     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
%     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
%     v3(n)=ftdcurve(p3,r3(n),1);
%     
%     %link 4 [s+1]
%     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
%     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
%     v4(n)=ftdcurve(p4,r4(n),1);
%     
%     %link 5 [s+2]
%     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
%     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
%     v5(n)=ftdcurve(p5,r5(n),1);
%     
%     if  (n>=horizon(1)) & (n<horizon(end))      
%         if (i<time+1)&(n==horizon(i))    
%             i=i+1
%             
%             %time delays
%             r1_1d = lag(r1(3:n),1,0);
%             r2_1d = lag(r2(3:n),1,0);
%             r3_1d = lag(r3(3:n),1,0);
%             r4_1d = lag(r4(3:n),1,0);
%             r5_1d = lag(r5(3:n),1,0);
%             q_ramp_1d=lag(q_ramp(3:n),1,0);
%             
%             %Define an F matrix as follows:
%             Z1=[r2(3:n) r1_1d r3_1d];            %for s-1 
%             Z2=[r3(3:n) r2_1d r4_1d q_ramp_1d];  %for s
%             Z3=[r4(3:n) r3_1d ];                 %for s+1 
%             
%             dtc=6; 
%             [Z1, m1, s1]=prepz(Z1, [], 1, 10, 0, dtc);  
%             [Z2, m2, s2]=prepz(Z2, [], 1, 10, 0, dtc);  
%             [Z3, m3, s3]=prepz(Z3, [], 1, 10, 0, dtc);  
%             z1=Z1; z2=Z2; z3=Z3;
%             
%             %             %============== link s-1
%             %             nn=[1 1 1 1 1 0]; [th, stats,e,var,Ps]=riv(z1, nn,[3 1 2 1],[0.71 0.25 0.04],[10000 10000 10000]); RT2=stats(3); YIC=stats(2); [at, bt]=getpar(th); 
%             %             %============== link s
%             %             nn= [1 1 1 1 1 1 1 0]; [th2, stats2,e2,var2,Ps2]=riv(z2, nn,[3 1 2 1],[0.76 0.15 0.12 0.15],[10000 10000 10000 10000 ]); RT2_2=stats2(3); YIC_2=stats2(2); [at_2, bt_2]=getpar(th2);
%             %             %=============== link s+1
%             %             nn=[1 1 1   0]; [th3, stats3,e3,var3,Ps3]=riv(z3, nn,[3 1 2 1],[0.78 0.24  ],[10000 10000 ]); RT2_3=stats3(3); YIC_3=stats3(2); [at_3, bt_3]=getpar(th3);
%             %============== link s-1     
%             Q1_up=cov(r2_1d); P1_us_l(i)=P1_up(i-1)+Q1_up; P1_up(i)=P1_us_l(i-1)-(P1_us_l(i-1)^2/(1+P1_us_l(i-1)));
%             Q2_up=cov(r1_1d); P2_us_l(i)=P2_up(i-1)+Q2_up; P2_up(i)=P2_us_l(i-1)-(P2_us_l(i-1)^2/(1+P2_us_l(i-1)));
%             Q3_up=cov(r3_1d); P3_us_l(i)=P3_up(i-1)+Q3_up; P3_up(i)=P3_us_l(i-1)-(P3_us_l(i-1)^2/(1+P3_us_l(i-1)));
%             nn=[1 1 1 1 1 0]; [th, stats,e,var,Ps]=riv(z1, nn,[1 0 0 0 0],[0.71 0.25 0.04],[P1_up(i) P2_up(i) P3_up(i)]); RT2=stats(3); YIC=stats(2); [at, bt]=getpar(th); 
%             
%             %============== link s
%             Q1_j=cov(r3_1d);      P1_j_l(i)=P1_j(i-1)+Q1_j; P1_j(i)=P1_j_l(i-1)-(P1_j_l(i-1)^2/(1+P1_j_l(i-1)));
%             Q2_j=cov(r2_1d);      P2_j_l(i)=P2_j(i-1)+Q2_j; P2_j(i)=P2_j_l(i-1)-(P2_j_l(i-1)^2/(1+P2_j_l(i-1)));
%             Q3_j=cov(r4_1d);      P3_j_l(i)=P3_j(i-1)+Q3_j; P3_j(i)=P3_j_l(i-1)-(P3_j_l(i-1)^2/(1+P3_j_l(i-1)));
%             Q4_j=cov(q_ramp_1d);  P4_j_l(i)=P4_j(i-1)+Q4_j; P4_j(i)=P4_j_l(i-1)-(P4_j_l(i-1)^2/(1+P4_j_l(i-1)));
%             nn= [1 1 1 1 1 1 1 0]; [th2, stats2,e2,var2,Ps2]=riv(z2, nn,[1 0 0 0 0],[0.76 0.15 0.12 0.15],[P1_j(i) P2_j(i) P3_j(i) P4_j(i)]); RT2_2=stats2(3); YIC_2=stats2(2); [at_2, bt_2]=getpar(th2);
%             
%             %=============== link s+1
%             Q1_dw=cov(r4_1d); P1_dw_l(i)=P1_dw(i-1)+Q1_dw; P1_dw(i)=P1_dw_l(i-1)-(P1_dw_l(i-1)^2/(1+P1_dw_l(i-1)));
%             Q2_dw=cov(r3_1d); P2_dw_l(i)=P2_dw(i-1)+Q2_dw; P2_dw(i)=P2_dw_l(i-1)-(P2_dw_l(i-1)^2/(1+P2_dw_l(i-1)));
%             nn=[1 1 1 0]; [th3, stats3,e3,var3,Ps3]=riv(z3, nn,[1 0 0 0 0],[0.78 0.24 ],[P1_dw(i) P2_dw(i)]); RT2_3=stats3(3); YIC_3=stats3(2); [at_3, bt_3]=getpar(th3);      
%             
%             %   ================ Extraction of model parameters =======================
%             aa(i,:)=[-at(2) -at_2(2) -at_3(2)];
%             bb(i,:)=[bt(1, 2) bt_2(1, 2) bt_3(1, 2)];
%             cc(i,:)=[0 bt_2(2, 2) 0];
%             dd(i,:)=[0 bt_2(3,2) 0];
%             
%             %Formulation of the NMSS matrix
%             F = [aa(i,1)  bb(i,1)      0  0;...
%                     bb(i,2)  aa(i,2)   cc(i,2) 0;...
%                     0     bb(i,3)   aa(i,3) 0;...
%                     -bb(i,2) -aa(i,2)  -cc(i,2) 1;];
%             g = [0 dd(i,2) 0 -dd(i,2)]'; 
%             %                    g = [0 0.15 0 -0.15]'; 
%             
%             
%             Q=eye(size(F)); Q(4,4)=10000000; Q(3,3)=1; Q(2,2)=1; Q(1,1)=1; 
%             r=1; k_adaptive(i,:) = dlqri(F,g,Q,r); 
%             
%             e(n) = uc(n) - r3(n); %self tuning controller
%             
%             q_ramp(n) = q_ramp(n-1) - k_adaptive(i,4)*e(n) - k_adaptive(i,1)*(r2(n)-r2(n-1)) - k_adaptive(i,2)*(r3(n)-r3(n-1)) - k_adaptive(i,3)*(r4(n)-r4(n-1));
%             
%             if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
%             if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
%             
%         else   
%             q_ramp(n)=q_ramp(n-1);
%         end
%     end  
%     if  n>=horizon(time)  
%         if N_ramp(n)>40; q_ramp(n)=16.7; %Release cars at capacity level
%         elseif (N_ramp(n)<=40)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
%         elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
%         end
%     end
%     waitbar(n/length(q1),h);
%     
% end
% close(h)
% 
% 
% figure(3); 
% subplot(211); h1=plot(r3,'r'); 
% subplot(212);plot(q_ramp,'r')
% % TRAFFIC Criteria
% TTT_adaptive = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(399)+q_ramp(399)-q4(399) ) * (960-400)/3600;
% TWT_adaptive = sum((q_urban(400:739)-q_ramp(400:739))) * (739-400)/3600;
% TTS_adaptive=TTT_adaptive+TWT_adaptive;
% 
% % % Control Criteria
% IAE_adaptive=sum(abs(r3(400:960)-uc(400:960)));
% resid=r3(400:960)-uc(400:960); index=find(resid>0); IAE2_adaptive=sum(resid(index));
% load_adaptive=sum(abs(q_ramp(400:960)));
% index=find(r3(400:960)<=uc(end)+0.05);Set_time_adaptive=index(1);
% 
% % figure(5);
% % subplot(221); plot(aa(:,2)); hold; plot(0.76*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('a(t)'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % subplot(222); plot(bb(:,2)); hold; plot(0.15*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('b(t)'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % subplot(223); plot(cc(:,2)); hold; plot(0.12*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('c(t)'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % subplot(224); plot(dd(:,2)); hold; plot(0.15*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('d(t)'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % 
% % figure(6);
% % subplot(221); plot(aa(:,1)); hold; plot(0.71*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('a(t) at s-1'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % subplot(222); plot(bb(:,1)); hold; plot(0.25*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('b(t) at s-1'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % subplot(223); plot(aa(:,3)); hold; plot(0.78*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('a(t) at s+1'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % subplot(224); plot(bb(:,3)); hold; plot(0.24*ones(94,1),'r');set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])
% % h1=ylabel('b(t) at s+1'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid
% % 
% % figure(4);
% % subplot(221); plot(k_adaptive(:,1));set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])  
% % h1=ylabel('Forward gain at s-1'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14); grid; hold; plot(1*ones(91,1),'r');
% % subplot(222); plot(k_adaptive(:,2));set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])  
% % h1=ylabel('Feedback gain at s'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14);grid; hold;plot(5*ones(91,1),'r');
% % subplot(223); plot(-k_adaptive(:,3));set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])  
% % h1=ylabel('Forward gain at s+1'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14);grid; hold;plot(0.8*ones(91,1),'r');
% % subplot(224); plot(k_adaptive(:,4));set(gca, 'xtick',[1 20 40 60 80 100], 'xticklabel', ['07:03'; '07:23'; '07:43'; '08:06'; '08:29'; '08:49'])  
% % h1=ylabel('Integral of error gain'); set(h1,'FontSize',14); h2=xlabel('Time (hr)'); set(h2,'FontSize',14);grid; hold;plot(-6.5*ones(91,1),'r');
% % 




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
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);%if (n>500)&(n<600); r5(n-1)=1.2*r5(n-1); end
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
    
    if  (n>=horizon(1)) & (n<horizon(end))      
        if (i<time+1)&(n==horizon(i))    
            i=i+1
            %81 times
            e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + 20*e(n);
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
            if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
        if N_ramp(n)>40; q_ramp(n)=16.7; %Release cars at capacity level
        elseif (N_ramp(n)<=40)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
        elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
        end
    end
    waitbar(n/length(q1),h);
    
end
close(h)

figure(3); subplot(211);h1=plot(r3,'m'); 
subplot(212);plot(q_ramp,'m');
% TRAFFIC Criteria
TTT_alinea = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(399)+q_ramp(399)-q4(399) ) * (960-400)/3600;
TWT_alinea = sum((q_urban(400:739)-q_ramp(400:739))) * (739-400)/3600;
TTS_alinea=TTT_alinea+TWT_alinea;

% % Control Criteria
IAE_alinea=sum(abs(r3(400:960)-uc(400:960)));
resid=r3(400:960)-uc(400:960); index=find(resid>0); IAE2_alinea=sum(resid(index));
load_alinea=sum(abs(q_ramp(400:960)));
index=find(r3(400:960)<=uc(end)+0.05);Set_time_alinea=index(1);


break

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SDP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1) PARAMETRIC
par_A=zeros(size(q1)); par_B=par_A; par_C=par_A;numMFs=[2:14];
e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); 

load sdp_results
% the following:
%state_A state_B state_C
%par_A par_B par_C 
% vx vx2 vx3 rules rules2 rules3 modsel modsel2 modsel3

h=waitbar(0,'Please wait.... Control in Progress (parametric) [5/7]');

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
            i=i+1
            e(n) = uc(n) - r3(n);
            
            % read-off States 
            state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x1_irw); state_A(n-2)=min(x1_irw); elseif state_A(n-2)>max(x1_irw); state_A(n-2)=max(x1_irw);end 
            state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x2_irw); state_B(n-2)=min(x2_irw); elseif state_B(n-2)>max(x2_irw); state_B(n-2)=max(x2_irw);end
            state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x3_irw); state_C(n-2)=min(x3_irw); elseif state_C(n-2)>max(x3_irw); state_C(n-2)=max(x3_irw);end
            
            % estimated parametric parameters 
            par_A(n-2)= FISnon_SISO(state_A(n-2), 1, numMFs(modsel),  [], VX,  rules,  {'gaussmf'});
            par_B(n-2)= FISnon_SISO(state_B(n-2), 1, numMFs(modsel2), [], VX2, rules2, {'gaussmf'});
            par_C(n-2)= FISnon_SISO(state_C(n-2), 1, numMFs(modsel3), [], VX3, rules3, {'gaussmf'});
            
            F = [0 0 0 0;
                par_B(n-2) par_A(n-2) par_C(n-2) 0;
                0 0 0 0;
                -par_B(n-2) -par_A(n-2) -par_C(n-2) 1];
            
            g =[0 pars(1,4) 0 -pars(1,4)];
            
            Q = eye(size(F)); Q(end,end)=10; r=1;  k=dlqri(F,g',Q,r); 
            
            k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
            
            e(n) = uc(n) - r3(n);  
            
            q_ramp(n)=q_ramp(n-1) + k_4(n)*e(n) - k_2(n)*(r3(n)-r3(n-1)) - k_1(n)*(r2(n)-r2(n-1)) - k_3(n)*(r4(n)-r4(n-1));  
            
            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
            
            if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
        if N_ramp(n)>40; q_ramp(n)=16.7; %Release cars at capacity level
        elseif (N_ramp(n)<=40)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
        elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
        end
    end
    waitbar(n/length(q1),h);
    
end

close(h)

figure(3); subplot(211);h1=plot(r3,'y'); legend('Open Loop','Critical','Reference Level','LLM-PIP','LLM-PIP2','ADAPTIVE','ALINEA','SDP');
subplot(212);plot(q_ramp,'y');
% TRAFFIC Criteria
TTT_sdp_par = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(399)+q_ramp(399)-q4(399) ) * (960-400)/3600;
TWT_sdp_par = sum((q_urban(400:739)-q_ramp(400:739))) * (739-400)/3600;
TTS_sdp_par=TTT_sdp_par+TWT_sdp_par;

% % Control Criteria
IAE_sdp_par=sum(abs(r3(400:960)-uc(400:960)));
resid=r3(400:960)-uc(400:960); index=find(resid>0); IAE2_sdp_par=sum(resid(index));
load_sdp_par=sum(abs(q_ramp(400:960)));
index=find(r3(400:960)<=uc(end)+0.05);Set_time_sdp_par=index(1);

% figure(4);
% subplot(311);hold; plot(state_A,  par_A,'go'); xlabel('Dependent State min(r(s,t),r(s-1,t-1))');   ylabel('Parameter a');
% subplot(312);hold; plot(state_B,  par_B,'go'); xlabel('Dependent State min(r(s-1,t),r(s,t-1))');   ylabel('Parameter b');
% subplot(313);hold; plot(state_C,  par_C,'go'); xlabel('Dependent State min(r(s+1,t),r(s+2,t-1))'); ylabel('Parameter c');


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SDP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %(2) Non-PARAMETRIC
% e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); 
% par_A_h=zeros(size(q1)); par_B_h=par_A_h; par_C_h=par_A_h; 
% 
% 
% h=waitbar(0,'Please wait.... Control in Progress (non-parametric) [6/7]');
% 
% for n=3:length(r2);
%     
%     i=1;
%     
%     %link 1 [s-2]
%     q1(n)=alpha1*min(v1(n-1),v2(n-1))*r1(n-1);
%     r1(n)=r1(n-1)+((q_in(n-1)-q1(n-1))*k1);
%     v1(n)=ftdcurve(p1,r1(n),1); 
%     
%     %link 2 [s-1]
%     q2(n)=alpha2*min(v2(n-1),v3(n-1))*r2(n-1);
%     r2(n)=r2(n-1)+((q1(n-1)-q2(n-1))*k2);
%     v2(n)=ftdcurve(p2,r2(n),1);
%     
%     %Ramp   
%     q_ramp(n)=alpha1*min(v_ramp(n-1),v3(n-1))*r_ramp(n-1);
%     r_ramp(n)=r_ramp(n-1)+((q_urban(n-1)-q_ramp(n-1))*kramp);    
%     v_ramp(n)=ftdcurve(p7,r_ramp(n),1);
%     N_ramp(n)=N_ramp(n-1)+alpha1*(q_urban(n-1)-q_ramp(n-1));
%     
%     
%     %link 3 [s]
%     q3(n)=alpha3*min(v3(n-1),v4(n-1))*r3(n-1);
%     r3(n)=r3(n-1)+((q_ramp(n-1)+q2(n-1)-q3(n-1))*k3);
%     v3(n)=ftdcurve(p3,r3(n),1);
%     
%     %link 4 [s+1]
%     q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);
%     r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4);
%     v4(n)=ftdcurve(p4,r4(n),1);
%     
%     %link 5 [s+2]
%     q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
%     r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
%     v5(n)=ftdcurve(p5,r5(n),1);
%     
%     if  (n>=horizon(1)) & (n<horizon(end))      
%         if (i<time+1)&(n==horizon(i))    
%             i=i+1
%             e(n) = uc(n) - r3(n);
%             
%             % read-off States 
%             state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x(:,1)); state_A(n-2)=min(x(:,1)); elseif state_A(n-2)>max(x(:,1)); state_A(n-2)=max(x(:,1));end 
%             state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x(:,2)); state_B(n-2)=min(x(:,2)); elseif state_B(n-2)>max(x(:,2)); state_B(n-2)=max(x(:,2));end
%             state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x(:,3)); state_C(n-2)=min(x(:,3)); elseif state_C(n-2)>max(x(:,3)); state_C(n-2)=max(x(:,3));end
%             
%             % estimated non-parametric parameters 
%             par_A_h(n-2) =  interp1(x1_irw,pars(:,1), state_A(n-2));
%             par_B_h(n-2) =  interp1(x2_irw,pars(:,2), state_B(n-2));
%             par_C_h(n-2) =  interp1(x3_irw,pars(:,3), state_C(n-2));
%             
%             F = [ 0 0 0 0;
%                 par_B_h(n-2) par_A_h(n-2) par_C_h(n-2) 0;
%                 0 0 0 0;
%                 -par_B_h(n-2) -par_A_h(n-2) -par_C_h(n-2) 1];
%             
%             g =[0 pars(1,4) 0 -pars(1,4)];
%             
%             Q = eye(size(F)); Q(end,end)=1; r=0.01;  k=dlqri(F,g',Q,r); 
%             
%             k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
%             
%             e(n) = uc(n) - r3(n);  
%             
%             q_ramp(n)=q_ramp(n-1) + k_4(n)*e(n) - k_2(n)*(r3(n)-r3(n-1)) - k_1(n)*(r2(n)-r2(n-1)) - k_3(n)*(r4(n)-r4(n-1));  
%             
%             if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>16.7; q_ramp(n)=16.7; end  
%             
%             if (n<450)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
%             
%         else   
%             q_ramp(n)=q_ramp(n-1);
%         end
%     end  
%     if  n>=horizon(time)  
%         if N_ramp(n)>40; q_ramp(n)=16.7; %Release cars at capacity level
%         elseif (N_ramp(n)<=40)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>16.7; q_ramp(n)=16.7; end   
%         elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>16.7; q_ramp(n)=16.7; end 
%         end
%     end
%     waitbar(n/length(q1),h);
%     
% end
% close(h)
% 
% figure(3);subplot(211); h1=plot(r3,'c'); 
% subplot(212);plot(q_ramp,'c')
% % TRAFFIC Criteria
% TTT_sdp_non_par = (sum(q_ramp(400:960)+q2(400:960)-q4(400:960)) + q2(399)+q_ramp(399)-q4(399) ) * (960-400)/3600;
% TWT_sdp_non_par = sum((q_urban(400:739)-q_ramp(400:739))) * (739-400)/3600;
% TTS_sdp_non_par=TTT_sdp_non_par+TWT_sdp_non_par;
% 
% % % Control Criteria
% IAE_sdp_non_par=sum(abs(r3(400:960)-uc(400:960)));
% resid=r3(400:960)-uc(400:960); index=find(resid>0); IAE2_sdp_non_par=sum(resid(index));
% load_sdp_non_par=sum(abs(q_ramp(400:960)));
% index=find(r3(400:960)<=uc(end)+0.05);Set_time_sdp_non_par=index(1);







ALINEA=['    [1]  ALINEA       '  num2str(IAE_alinea)  '     ' num2str(IAE2_alinea)  '     '  num2str(load_alinea) ...
        '       ' num2str(load_alinea) '               ' num2str(Set_time_alinea)];
LLM   =['    [2]  LLM          '  num2str(IAE_llm)  '      ' num2str(IAE2_llm)  '    '  num2str(load_llm) ...
        '       ' num2str(load_llm) '               ' num2str(Set_time_llm)];
LLM_d =['    [3]  LLM_det      '  num2str(IAE_llm2)  '      ' num2str(IAE2_llm2)  '    '  num2str(load_llm2) ...
        '       ' num2str(load_llm2) '                ' num2str(Set_time_llm2)];
ADAPT =['    [4]  ADAPTIVE     '  num2str(IAE_adaptive)  '      ' num2str(IAE2_adaptive)  '    '  num2str(load_adaptive) ...
        '       ' num2str(load_adaptive) '               ' num2str(Set_time_adaptive)];
SDPpar=['    [5]  SDP-par      '  num2str(IAE_sdp_par)  '      ' num2str(IAE2_sdp_par)  '    '  num2str(load_sdp_par)  ...
        '       ' num2str(load_sdp_par) '               ' num2str(Set_time_sdp_par)];
% SDP   =['    [6]  SDP-non.par  '  num2str(IAE_sdp)  '      ' num2str(IAE2_sdp)  '    '  num2str(load_sdp) ...
%         '        '  num2str(load_sdp) '               ' num2str(Set_time_sdp_npar)];

% SDPdet=['    [7]  SDP-Det      '  num2str(IAE_det)  '      ' num2str(IAE2_det)  '    '  num2str(load_det) ...
%         '       ' num2str(load_det) '               ' num2str(Set_time_sdp_det)];


fprintf('                   Control Criteria                    Results                      \n \n');
fprintf('      Controller       IAE        Sum(r +)   Ramp Effort     Load Reduct       Settling time    \n  ');
disp('---------------------------------------------------------------------------------------------------');
disp(ALINEA)
disp(LLM)
disp(LLM_d)
disp(ADAPT)
% disp(SDP)
disp(SDPpar)
% disp(SDPdet)




ALINEA=['    [1]  ALINEA       '  num2str(TTT_alinea)  '      ' num2str(TWT_alinea)  '    '  num2str(TTS_alinea)];
LLM   =['    [2]  LLM          '  num2str(TTT_llm)  '      ' num2str(TWT_llm)  '    '  num2str(TTS_llm) ];
LLM_d =['    [3]  LLM_det      '  num2str(TTT_llm2)  '      ' num2str(TWT_llm2)  '    '  num2str(TTS_llm2)];
ADAPT =['    [4]  ADAPTIVE     '  num2str(TTT_adaptive)  '      ' num2str(TWT_adaptive)  '    '  num2str(TTS_adaptive) ];
SDPpar=['    [5]  SDP-par      '  num2str(TTT_sdp_par)  '      ' num2str(TWT_sdp_par)  '    '  num2str(TTS_sdp_par) ];
% SDP   =['    [6]  SDP-non.par  '  num2str(IAE_sdp)  '      ' num2str(IAE2_sdp)  '    '  num2str(load_sdp) ...
%         '        '  num2str(load_sdp) '               ' num2str(Set_time_sdp_npar)];

% SDPdet=['    [7]  SDP-Det      '  num2str(IAE_det)  '      ' num2str(IAE2_det)  '    '  num2str(load_det) ...
%         '       ' num2str(load_det) '               ' num2str(Set_time_sdp_det)];


fprintf('                  Traffic Criteria:                    Results                      \n \n');
fprintf('      Controller       TTS           TWT          TTS    \n  ');
disp('---------------------------------------------------------------------------------------------------');
disp(ALINEA)
disp(LLM)
disp(LLM_d)
disp(ADAPT)
% disp(SDP)
disp(SDPpar)
% disp(SDPdet)



