close all; clear all; clc
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

% initiali4se other variables
r1=NaN*ones(size(q_in)); r2=NaN*ones(size(q_in)); r3=NaN*ones(size(q_in)); r4=NaN*ones(size(q_in)); r5=NaN*ones(size(q_in));  
v1=NaN*ones(size(q_in)); v2=NaN*ones(size(q_in)); v3=NaN*ones(size(q_in)); v4=NaN*ones(size(q_in)); v5=NaN*ones(size(q_in)); 
q1=NaN*ones(size(q_in)); q2=NaN*ones(size(q_in)); q3=NaN*ones(size(q_in)); q4=NaN*ones(size(q_in)); q5=NaN*ones(size(q_in)); q_ramp=NaN*ones(size(q_in));
r1(1:2)=1; r2(1:2)=2; r3(1:2)=1; r4(1:2)=1; r5(1:2)=1; q1(1:2)=1; q2(1:2)=1; q3(1:2)=1; q4(1:2)=1; q5(1:2)=1; 
v1(1:2)=ftdcurve(p1,r1(1:2),1); v2(1:2)=ftdcurve(p2,r2(1:2),1); v3(1:2)=ftdcurve(p3,r3(1:2),1); v4(1:2)=ftdcurve(p4,r4(1:2),1);
v5(1:2)=ftdcurve(p5,r5(1:2),1); dt=1; r_ramp(1:2)=1; v_ramp(1:2)=ftdcurve(p7,r_ramp(1:2),1);; q_ramp(1:2)=1;
r3_llm=NaN*ones(size(r2));r3_llm(1:2)=r3(1:2);
N_ramp=zeros(size(q_in)); N_ramp(1:3)=q_ramp(1:3);

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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1);if (n>1620)&(n<1680); q4(n)=1/3*q4(n); end
    r4(n)=r4(n-1)+((q3(n-1)-q4(n-1))*k4); 
    v4(n)=ftdcurve(p4,r4(n),1);
    
    %link 5 [s+2]
    q5(n)=alpha5*min(v5(n-1),v6(n-1))*r5(n-1);
    r5(n)=r5(n-1)+((q4(n-1)-q5(n-1))*k5);    
    v5(n)=ftdcurve(p5,r5(n),1);

     r3_llm(n)=bb(2)*r2(n-1)+aa(2)*r3(n-1)+cc(2)*r4(n-1)+dd(2)*q_ramp(n-1);
     
    waitbar(n/length(r2),h);
end

close(h)

figure(3);
subplot(211);plot(r3,'k'); hold;plot(p3(4)*ones(length(q1),1),'-.'); 
h3=ylabel('Density (veh/km)');           set(h3,'FontSize',14);
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'; '12:00'])  
% subplot(211);h1=plot(r4,'k'); hold;plot(p4(4)*ones(length(q1),1),'-.'); 
h3=ylabel('Density (veh/km)');           set(h3,'FontSize',14);
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
%set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'; '12:00'])  
subplot(212);hold; plot(q_urban,'k'); plot(q_ramp,'k');
h3=ylabel(' Flow (veh/min)'); set(h3,'FontSize',14);
h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
%set(gca, 'xtick',[1 360 720 1080 1440 1800 2160], 'xticklabel', ['06:00'; '07:00'; '08:00'; '09:00'; '10:00'; '11:00'; '12:00'])  

% figure(11); 
% subplot(221);plot(q2);hold
% subplot(222);plot(q3);hold
% subplot(223);plot(q4);hold
% subplot(224);plot(q5);hold

%time delays
r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r1)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4)); r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5));
v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v4)); v4_2d = lag(v4,2,mean(v4)); v5_2d = lag(v5,2,mean(v5));
q_ramp_1d=lag(q_ramp,1,mean(q_ramp));

figure(6);
subplot(321); plot(r1,q1,'.'); hold; plot(p1(4),linspace(0,120),'k.'); plot(p1(4),max(q1),'r.');title('s-2');
subplot(322); plot(r2,q2,'.'); hold; plot(p2(4),linspace(0,120),'k.'); plot(p2(4),max(q2),'r.');title('s-1 {load}');
subplot(323); plot(r3,q3,'.'); hold; plot(p3(4),linspace(0,120),'k.'); plot(p3(4),max(q3),'r.');title('s {junction}');
subplot(324); plot(r4,q4,'.'); hold; plot(p4(4),linspace(0,120),'k.'); plot(p4(4),max(q4),'r.');title('s+1');
subplot(325); plot(r5,q5,'.'); hold; plot(p5(4),linspace(0,120),'k.'); plot(p5(4),max(q5),'r.');title('s+2');
subplot(326); plot(r_ramp,q_ramp,'.'); hold; plot(p7(4),linspace(0,120),'k.'); plot(p7(4),max(q_ramp),'r.');title('ramp');
xlabel('Density veh/km'); ylabel('Flow veh/h');

figure(4);
x=[0:1:500]';y_1=ftdcurve(p1,x,1);subplot(231);plot(x,y_1); hold; axis([0 200 0 130]);title('s-2');
x=[0:1:500]';y_2=ftdcurve(p2,x,1);subplot(232);plot(x,y_2); hold; axis([0 200 0 130]);title('s-1 {load}');
x=[0:1:500]';y_3=ftdcurve(p3,x,1);subplot(233);plot(x,y_3); hold; axis([0 200 0 130]);title('s {junction}');
x=[0:1:500]';y_4=ftdcurve(p4,x,1);subplot(234);plot(x,y_4); hold; axis([0 200 0 130]);title('s+1'); 
x=[0:1:500]';y_5=ftdcurve(p5,x,1);subplot(235);plot(x,y_5); hold; axis([0 200 0 130]);title('s+2'); 
x=[0:1:500]';y_6=ftdcurve(p7,x,1);subplot(236);plot(x,y_6); hold; axis([0 375 0 100]);title('ramp');
subplot(231);plot(r1,v1,'.'); plot(p1(4),p1(2),'c.');title('s-2');xlabel('Density');ylabel('Velocity');
subplot(232);plot(r2,v2,'.'); plot(p2(4),p2(2),'c.');title('s-1');xlabel('Density');ylabel('Velocity');
subplot(233);plot(r3,v3,'.'); plot(p3(4),p3(2),'c.');title('s');  xlabel('Density');ylabel('Velocity');
subplot(234);plot(r4,v4,'.'); plot(p4(4),p4(2),'c.');title('s+1');xlabel('Density');ylabel('Velocity');
subplot(235);plot(r5,v5,'.'); plot(p5(4),p5(2),'c.');title('s+2');xlabel('Density');ylabel('Velocity'); 
subplot(236);plot(r_ramp,v_ramp,'.');plot(p7(4),p7(2),'c.');title('ramp');xlabel('Density');ylabel('Velocity');
 
figure(5); 
o1=occ2den(r1,4,2,3,0);     subplot(321); plot(q1,o1,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o2=occ2den(r2,4,2,3,0);     subplot(322); plot(q2,o2,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o3=occ2den(r3,4,2,3,0);     subplot(323); plot(q3,o3,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o4=occ2den(r4,4,2,3,0);     subplot(324); plot(q4,o4,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o5=occ2den(r5,4,2,3,0);     subplot(325); plot(q5,o5,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o6=occ2den(r_ramp,4,2,3,0); subplot(326); plot(q_ramp,o6,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');hold;

% -------------- Section 3 Conrol Implementation -------------
uc=[zeros(600,1); (p3(4)+5)*ones(300,1); (p3(4)-5)*ones(300,1); (p3(4)+8)*ones(300,1); p3(4)*ones(500,1)];
start=500; finish=2000;figure(3); subplot(211); plot(uc,'r'); horizon=start:6:finish; time=length(horizon);
%------------------------------------------------------------------

% *** LLM Controller on STM model ***
load llm_par

%Formulation of the NMSS matrix
F = [aa(1)  bb(1)      0   0;...
     bb(2)  aa(2)   cc(2)  0;...
        0   bb(3)   aa(3)  0;...
    -bb(2) -aa(2)  -cc(2)  1;];

g = [0 dd(2) 0 -dd(2)]'; 

N_ramp=zeros(size(q_in)); N_ramp(1:3)=q_ramp(1:3);

% F = [aa(1)  bb(1)      0   0;...
%      bb(2)  aa(2)   cc(2)  0;...
%         0   bb(3)   aa(3)  0;...
%         0  -bb(3)  -aa(3)  1;];
% 
% g = [0 dd(2) 0 0]'; 

Q=eye(size(F)); Q(4,4)=10000; Q(3,3)=1; Q(2,2)=1; Q(1,1)=1; 

r=1; k_llm = dlqri(F,g,Q,r); k_llm(end)=-k_llm(end); 
max(r2)/p2(4)

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
    q4(n)=alpha4*min(v4(n-1),v5(n-1))*r4(n-1); if (n>1660)&(n<1720); q4(n)=1/3*q4(n); end
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
            
            if (n<900)&(q_ramp(n)>q_urban(n)); q_ramp(n)=q_urban(n); end
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)  
        if N_ramp(n)>100; q_ramp(n)=20; %Release cars at capacity level
        elseif (N_ramp(n)<=100)&(N_ramp(n)>=0); q_ramp(n)=0.5*N_ramp(n); if q_ramp(n)>20; q_ramp(n)=20; end   
        elseif N_ramp(n)<10; q_ramp(n)=q_urban(n);   if q_ramp(n)>20; q_ramp(n)=20; end 
        end
    end
    waitbar(n/length(q1),h);
    
end
close(h)

figure(3); 
subplot(211);h1=plot(r3,'g'); 
subplot(212); h1=plot(q_ramp,''); %h5=legend('Urban Demand', 'Unmetered Ramp Flow', 'Metered Ramp Flow');


figure(5); 
o1=occ2den(r1,4,2,3,0);     subplot(321); plot(q1,o1,'r.');     %ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o2=occ2den(r2,4,2,3,0);     subplot(322); plot(q2,o2,'r.');     %ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o3=occ2den(r3,4,2,3,0);     subplot(323); plot(q3,o3,'r.');     %ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o4=occ2den(r4,4,2,3,0);     subplot(324); plot(q4,o4,'r.');     %ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o5=occ2den(r5,4,2,3,0);     subplot(325); plot(q5,o5,'r.');     %ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold;
o6=occ2den(r_ramp,4,2,3,0); subplot(326); plot(q_ramp,o6,'r.'); %ylabel('Occupancy %'); xlabel('Flow veh/sample');hold;
% 

figure(6);
subplot(321); plot(r1,q1,'r.'); 
subplot(322); plot(r2,q2,'r.'); 
subplot(323); plot(r3,q3,'r.');
subplot(324); plot(r4,q4,'r.');
subplot(325); plot(r5,q5,'r.');
subplot(326); plot(r_ramp,q_ramp,'.');
