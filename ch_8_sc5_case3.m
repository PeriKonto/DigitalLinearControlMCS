% case 3 scenario 5 
% this addresses the question adaptive control on SDP-PIP as compared to the normal fixed gain controlclear all; close all; clc
clear al; close all;
load llm_par

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
p7=[2   50   60  25   3600];   %ramp

%FACTORS AFFECTING TRAFFIC CONGESTION SEVERITY
X=0.6; Y=1; Z=1.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%% boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t, q, v, o]=m27load(26, 10, 1, 1); r=occ2den(o, 4, 2, 3, 1); 
q_in=[X*q(:,1)]; q_in=irwsm(q_in,1,0.001); q_in=resample(q_in,6,1); 
v6=irwsm(v(:,7),1,0.001); v6=resample(v6,6,1); 
v6=[110*ones(100,1); Y*v6(2000:4250)]; q_in=[10*ones(100,1); q_in(2160:4320)]; 
q_urban=[q(:,2)]; q_urban=irwsm(q_urban,1,0.001);  q_urban=Z*resample(q_urban,6,1);q_urban=[2*ones(100,1); q_urban(2160:4320)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise other variables
r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
r3_llm=NaN*ones(size(r2));r3_llm(1:2)=r3(1:2); N_ramp=zeros(size(q_in)); N_ramp(1:3)=q_ramp(1:3);

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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1); if ((n>900)&(n<1000)); q4(n)=1/3*q4(n); end 
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4); 
    v4(n)=ftdcurve(p4,r4(n),1);
    
    %link 5 [s+2]
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);
     
    waitbar(n/length(r2),h);
end

close(h)

figure(3);subplot(211);plot(r3,''); hold; plot(p3(4)*ones(length(q1),1),'-.'); h3=ylabel('Density (veh/km)'); set(h3,'FontSize',14); h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
% set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'; '12:00'])  
subplot(212);hold; plot(q_urban,'k'); plot(q_ramp,'k'); h3=ylabel(' Flow (veh/min)'); set(h3,'FontSize',14);
h4=xlabel('Time (hr)'); set(h4,'FontSize',14); 
% set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'; '12:00'])  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% time delays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r1)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4)); r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5));
v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v4)); v4_2d = lag(v4,2,mean(v4)); v5_2d = lag(v5,2,mean(v5));
q_ramp_1d=lag(q_ramp,1,mean(q_ramp));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(6);
% subplot(321); plot(r1,q1,'.'); hold; plot(p1(4),linspace(0,120),'k.'); plot(p1(4),max(q1),'r.');title('s-2');
% subplot(322); plot(r2,q2,'.'); hold; plot(p2(4),linspace(0,120),'k.'); plot(p2(4),max(q2),'r.');title('s-1 {load}');
% subplot(323); plot(r3,q3,'.'); hold; plot(p3(4),linspace(0,120),'k.'); plot(p3(4),max(q3),'r.');title('s {junction}');
% subplot(324); plot(r4,q4,'.'); hold; plot(p4(4),linspace(0,120),'k.'); plot(p4(4),max(q4),'r.');title('s+1');
% subplot(325); plot(r5,q5,'.'); hold; plot(p5(4),linspace(0,120),'k.'); plot(p5(4),max(q5),'r.');title('s+2');
% subplot(326); plot(r_ramp,q_ramp,'.'); hold; plot(p7(4),linspace(0,120),'k.'); plot(p7(4),max(q_ramp),'r.');title('ramp');
% xlabel('Density veh/km'); ylabel('Flow veh/h');
% 
% figure(5); 
% o1=occ2den(r1,4,2,3,0);     subplot(321); plot(q1,o1,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q1),linspace(0,20),'r')
% o2=occ2den(r2,4,2,3,0);     subplot(322); plot(q2,o2,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q2),linspace(0,20),'r')
% o3=occ2den(r3,4,2,3,0);     subplot(323); plot(q3,o3,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q3),linspace(0,20),'r')
% o4=occ2den(r4,4,2,3,0);     subplot(324); plot(q4,o4,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q4),linspace(0,20),'r')
% o5=occ2den(r5,4,2,3,0);     subplot(325); plot(q5,o5,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q5),linspace(0,20),'r')
% o6=occ2den(r_ramp,4,2,3,0); subplot(326); plot(q_ramp,o6,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');hold;
% traffic criteria ... 

index2=find(r3==max(r3));
index=find(r3>p3(4));

%Period density is above critical
TTT_open = (sum(q_ramp(index(1):index2)+q2(index(1):index2)-q5(index(1):index2)) +...
            q2(index(1)-1)+q_ramp(index(1)-1)-q5(index(1)-1) ) * (index2-index(1))/3600;

value2=index(1);  %activation of control
finish=1200; %deactivation of control 

uc=[zeros(500,1); (p3(4)+2)*ones(1800,1)];

horizon=value2:6:finish; time=length(horizon);


% *** LLM Controller on STM model ***
load llm_par

%Formulation of the NMSS matrix
F = [aa(1)  bb(1)      0   0;...
     bb(2)  aa(2)   cc(2)  0;...
        0   bb(3)   aa(3)  0;...
    -bb(2) -aa(2)  -cc(2)  1;];
g = [0 dd(2) 0 -dd(2)]'; 

N_ramp=zeros(size(q_in)); N_ramp(1:3)=q_ramp(1:3); Q=eye(size(F)); 

c1=1; c2=1; c3=100;  c4=300;

Qb=1; C=[0 0 c3 c4]; Q=C'*Qb*C; %Q=Q+eye(4); 

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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);if ((n>900)&(n<1000)); q4(n)=1/3*q4(n); end 
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
        
         q_ramp(n)=min(((p3(4)-15)/r3(n))*20,20);   if q_ramp(n)>20; q_ramp(n)=20; elseif N_ramp(n)<15; q_ramp(n)=q_urban(n); end   
    end

    waitbar(n/length(q1),h);
    
end
close(h)

figure(3); subplot(211);h1=plot(r3,'k');plot(uc,'r'); subplot(212); h1=plot(q_ramp,'k'); 

IAE_llm=sum(abs(r3(value2:finish)-56)); 
resid_1=r3(value2:finish)-56; count=find(resid_1>0); IAE2_llm=sum(resid_1(count));
TTT_llm=(sum(q_ramp(value2:1600)+q2(value2:1600)-q5(value2:1600)) + q2(value2-1) ...
                       +q_ramp(value2-1)-q5(value2-1) ) * (1600-value2)/3600;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par_A=zeros(size(q1)); par_B=par_A; par_C=par_A;numMFs=[2:14]; e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); 
error=zeros(size(q1));

load sdp_results

h=waitbar(0,'Please wait.... Control in Progress (parametric) [5/7]');
 c4_n=zeros(size(time)); c3_n=c4_n;
 
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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1); if ((n>900)&(n<1000)); q4(n)=1/3*q4(n); end 
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
            
            % read-off States 
            state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x1_irw); state_A(n-2)=min(x1_irw); elseif state_A(n-2)>max(x1_irw); state_A(n-2)=max(x1_irw);end 
            state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x2_irw); state_B(n-2)=min(x2_irw); elseif state_B(n-2)>max(x2_irw); state_B(n-2)=max(x2_irw);end
            state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x3_irw); state_C(n-2)=min(x3_irw); elseif state_C(n-2)>max(x3_irw); state_C(n-2)=max(x3_irw);end
            
            % estimated parametric parameters 
            par_A(n-2)= FISnon_SISO(state_A(n-2), 1, numMFs(modsel),  [], VX,  rules,  {'gaussmf'});
            par_B(n-2)= FISnon_SISO(state_B(n-2), 1, numMFs(modsel2), [], VX2, rules2, {'gaussmf'});
            par_C(n-2)= FISnon_SISO(state_C(n-2), 1, numMFs(modsel3), [], VX3, rules3, {'gaussmf'});
            
            F = [1 0 0 0;
                par_B(n-2) par_A(n-2) par_C(n-2) 0;
                0 0 1 0;
                -par_B(n-2) -par_A(n-2) -par_C(n-2) 1];
            
            g =[0 pars(1,4) 0 -pars(1,4)];
            
            c3_n(i)=(r4(n)/p3(4))*(r4(n)/10); 
                                                            
            c1=1; c2=1; c3=c3_n(i); c4=10;%c4_n(i);
            
            Qb=1; C=[c1 c2 c3 c4]; Q=C'*Qb*C; Q=Q+eye(4); 
            
            r=1; k=dlqri(F,g',Q,r);
            
            k_1(i)=k(1,1); k_2(i) = k(1,2);  k_3(i) = k(1,3);  k_4(i) = -k(1,4); 
            
            e(n) = uc(n) - r3(n);
            
            q_ramp(n)=q_ramp(n-1) + k_4(i)*e(n) - k_2(i)*(r3(n)-r3(n-1)) - k_1(i)*(r2(n)-r2(n-1)) - k_3(i)*(r4(n)-r4(n-1));  
            
           if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)
        
         q_ramp(n)=min(((p3(4)-15)/r3(n))*20,20);   if q_ramp(n)>20; q_ramp(n)=20; elseif N_ramp(n)<15; q_ramp(n)=q_urban(n); end   
         
    end

    waitbar(n/length(q1),h);
    
end

close(h)

figure(3); subplot(211);h1=plot(r3,'m'); subplot(212);plot(q_ramp,'m');

IAE_sdp=sum(abs(r3(value2:finish)-56)); 
resid_1=r3(value2:1600)-56; count=find(resid_1>0); IAE2_sdp=sum(resid_1(count));
TTT_sdp=(sum(q_ramp(value2:1600)+q2(value2:1600)-q5(value2:1600)) + q2(value2-1) ...
                       +q_ramp(value2-1)-q5(value2-1) ) * (1600-value2)/3600;
