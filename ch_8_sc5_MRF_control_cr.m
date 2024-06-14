% scenario 5
% Demonstration of SDP-PIP properties compared to a fixed gain control.
%multiple reference levels
clear all; close all; clc

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
X=0.67; Y=1; Z=1.8;

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

figure(3);
subplot(311);plot(r3,'k'); hold; plot(p3(4)*ones(length(q1),1),'-.'); 
h3=ylabel('Density (veh/km)'); set(h3,'FontSize',14); h4=xlabel('Time (hr)'); set(h4,'FontSize',14);
subplot(312);plot(r3,'k'); hold; plot(p3(4)*ones(length(q1),1),'-.'); 
subplot(313);hold; plot(q_urban,'k'); plot(q_ramp,'k'); h3=ylabel(' Flow (veh/min)'); set(h3,'FontSize',14);
h4=xlabel('Time (hr)'); set(h4,'FontSize',14); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% time delays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1_1d = lag(r1,1,mean(r1)); r1_2d = lag(r1,2,mean(r1)); r2_1d = lag(r2,1,mean(r1)); r2_2d = lag(r2,2,mean(r2)); r3_1d = lag(r3,1,mean(r3)); r3_2d = lag(r3,2,mean(r3)); r4_1d = lag(r4,1,mean(r4)); r4_2d = lag(r4,2,mean(r4)); r5_1d = lag(r5,1,mean(r5)); r5_2d = lag(r5,2,mean(r5));
v2_2d = lag(v2,2,mean(v2)); v3_1d = lag(v3,1,mean(v3)); v3_2d = lag(v3,2,mean(v3)); v4_1d = lag(v4,1,mean(v4)); v4_2d = lag(v4,2,mean(v4)); v5_2d = lag(v5,2,mean(v5));
q_ramp_1d=lag(q_ramp,1,mean(q_ramp));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);
subplot(321); plot(r1,q1,'.'); hold; plot(p1(4),linspace(0,120),'k.'); plot(p1(4),max(q1),'r.');title('s-2');
subplot(322); plot(r2,q2,'.'); hold; plot(p2(4),linspace(0,120),'k.'); plot(p2(4),max(q2),'r.');title('s-1 {load}');
subplot(323); plot(r3,q3,'.'); hold; plot(p3(4),linspace(0,120),'k.'); plot(p3(4),max(q3),'r.');title('s {junction}');
subplot(324); plot(r4,q4,'.'); hold; plot(p4(4),linspace(0,120),'k.'); plot(p4(4),max(q4),'r.');title('s+1');
subplot(325); plot(r5,q5,'.'); hold; plot(p5(4),linspace(0,120),'k.'); plot(p5(4),max(q5),'r.');title('s+2');
subplot(326); plot(r_ramp,q_ramp,'.'); hold; plot(p7(4),linspace(0,120),'k.'); plot(p7(4),max(q_ramp),'r.');title('ramp');
xlabel('Density veh/km'); ylabel('Flow veh/h');

figure(5); 
o1=occ2den(r1,4,2,3,0);     subplot(321); plot(q1,o1,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q1),linspace(0,20),'r')
o2=occ2den(r2,4,2,3,0);     subplot(322); plot(q2,o2,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q2),linspace(0,20),'r')
o3=occ2den(r3,4,2,3,0);     subplot(323); plot(q3,o3,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q3),linspace(0,20),'r')
o4=occ2den(r4,4,2,3,0);     subplot(324); plot(q4,o4,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q4),linspace(0,20),'r')
o5=occ2den(r5,4,2,3,0);     subplot(325); plot(q5,o5,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');    hold; plot(max(q5),linspace(0,20),'r')
o6=occ2den(r_ramp,4,2,3,0); subplot(326); plot(q_ramp,o6,'.'); ylabel('Occupancy %'); xlabel('Flow veh/sample');hold;

imposed_load=max(r2)/p2(4);
imposed_load2=max(r5)/p5(4);
index=find(r3>p3(4));

% -------------- Section 3 Conrol Implementation -------------
uc=[zeros(600,1); (p3(4)-5)*ones(300,1); (p3(4)+5)*ones(300,1); (p3(4)+10)*ones(200,1); (p3(4)+2)*ones(800,1)];
start=600; finish=2100;figure(3); 
subplot(311); plot(uc,'r'); 
subplot(312); plot(uc,'r'); 

horizon=start:6:finish; time=length(horizon);
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

Q=eye(size(F)); Q(4,4)=10000; Q(3,3)=1; Q(2,2)=1; Q(1,1)=1; 

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
        
         q_ramp(n)=min((p3(4)/r3(n))*10,20);   if q_ramp(n)>20; q_ramp(n)=20; end   
    end

    waitbar(n/length(q1),h);
    
end
close(h)

figure(3); 
subplot(311); plot(r3(1:2100),'g'); 
subplot(312); plot(r3(1:2100),'g');
subplot(313); plot(q_ramp(1:2100),'g'); 

figure(6); subplot(321); plot(r1,q1,'r.'); subplot(322); plot(r2,q2,'r.'); subplot(323); plot(r3,q3,'r.'); subplot(324); plot(r4,q4,'r.'); subplot(325); plot(r5,q5,'r.'); subplot(326); plot(r_ramp,q_ramp,'.');

figure(5); 
o1=occ2den(r1,4,2,3,0); subplot(321); plot(q1,o1,'r.'); o2=occ2den(r2,4,2,3,0);     subplot(322); plot(q2,o2,'r.');
o3=occ2den(r3,4,2,3,0); subplot(323); plot(q3,o3,'r.'); o4=occ2den(r4,4,2,3,0);     subplot(324); plot(q4,o4,'r.');
o5=occ2den(r5,4,2,3,0); subplot(325); plot(q5,o5,'r.'); o6=occ2den(r_ramp,4,2,3,0); subplot(326); plot(q_ramp,o6,'r.');

imposed_load_llm=max(r2)/p2(4);

IAE_llm_1=sum(abs(r3(start:start+300)-uc(start:start+300))); 
IAE_llm_2=sum(abs(r3(start+300:start+600)-uc(start+300:start+600))); 
IAE_llm_3=sum(abs(r3(start+600:start+800)-uc(start+600:start+800))); 
IAE_llm_4=sum(abs(r3(start+800:1630)-uc(start+800:1630))); 

resid_1=r3(start:start+300)-uc(start:start+300);         count=find(resid_1>0);   IAE2_llm_1=sum(resid_1(count));
resid_2=r3(start+300:start+600)-uc(start+300:start+600); count2=find(resid_2>0);  IAE2_llm_2=sum(resid_2(count2));
resid_3=r3(start+600:start+800)-uc(start+600:start+800); count3=find(resid_3>0);  IAE2_llm_3=sum(resid_3(count3));
resid_4=r3(start+800:1630)-uc(start+800:1630);           count4=find(resid_4>0);  IAE2_llm_4=sum(resid_4(count4));

load_llm_1=sum(q_ramp(start:start+300))/6;
load_llm_2=sum(q_ramp(start+300:start+600))/6;
load_llm_3=sum(q_ramp(start+600:start+800))/6;
load_llm_4=sum(q_ramp(start+800:1630))/6;

indexs =find( (r3(start:start+300)<=uc(start+300)+0.01)     & (r3(start:start+300)>uc(start+300)-0.01) );    count=indexs-lag(indexs,1,0);    I=find(count==1);   Set_time_llm_1=indexs(I(1));  
indexs2=find( (r3(start+300:start+600)<=uc(start+600)+0.01) & (r3(start+300:start+600)>uc(start+600)-0.01)); count2=indexs2-lag(indexs2,1,0); I2=find(count2==1); Set_time_llm_2=indexs2(I2(1));  
indexs3=find( (r3(start+600:start+800)<=uc(start+800)+0.01) & (r3(start+600:start+800)>uc(start+800)-0.01)); count3=indexs3-lag(indexs3,1,0); I3=find(count3==1); Set_time_llm_3=indexs3(I3(1));
indexs4=find( (r3(start+800:1630)<=uc(1630)+0.01)           & (r3(start+800:1630)<=uc(1630)-0.01));          count4=indexs4-lag(indexs4,1,0); I4=find(count4==1); Set_time_llm_4=indexs4(I4(1)); 



% % *** LLM Controller on STM model ***
% load llm2_par
% % Formulation of the NMSS matrices
% F_llm2=[ 1   aa2(1)  1  cc2(1)  0     0       0 ; ...
%          1    0      0   0      0     0       0 ; ...     
%          0   bb2(2)  1  aa2(2)  1   cc2(2)    0 ; ...
%          0     0     1   0      0     0       0 ; ...
%          0     0     0  bb2(3)  1   aa2(3)    0 ; ...   
%          0     0     0   0      1     0       0 ; ...
%          0  -bb2(2) -1 -aa2(2) -1  -cc2(2)    1 ];
% 
% g_llm2 =[0 0 dd2(1) 0 0 0 -dd2(1)]';
% 
% Q=eye(size(F_llm2)); Q(end,end)=1000;
% Q(1,1)=1; Q(2,2)=1; Q(3,3)=1; Q(4,4)=1; Q(5,5)=1; Q(6,6)=1;
% 
% r=1; k_llm2 = dlqri(F_llm2,g_llm2,Q,r); k_llm2(end)=-k_llm2(end); 
% 
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
%             i=i+1;
%           
%             e(n) = uc(n) - r3(n);
%             
%             q_ramp(n) = q_ramp(n-1) + k_llm2(7)*e(n) ...
%                                     - k_llm2(1)*(r2(n)-r2(n-1)) ...
%                                     - k_llm2(2)*(r2(n-1)-r2(n-2)) ...
%                                     - k_llm2(3)*(r3(n)-r3(n-1)) ...
%                                     - k_llm2(4)*(r3(n-1)-r3(n-2)) ...
%                                     - k_llm2(5)*(r4(n)-r4(n-1)) ...
%                                     - k_llm2(6)*(r4(n-1)-r4(n-2));  
%             
%            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
%             
%         else   
%             q_ramp(n)=q_ramp(n-1);
%         end
%     end  
%     if  n>=horizon(time)
%         
%          q_ramp(n)=min((p3(4)/r3(n))*10,20);   if q_ramp(n)>20; q_ramp(n)=20; end   
%     end
% 
%     waitbar(n/length(q1),h);
%     
% end
% close(h)
% 
% figure(3); subplot(211);h1=plot(r3,'k'); subplot(212);plot(q_ramp,'k');
% 
% IAE_llm2_1=sum(abs(r3(start:start+300)-uc(start:start+300))); 
% IAE_llm2_2=sum(abs(r3(start+300:start+600)-uc(start+300:start+600))); 
% IAE_llm2_3=sum(abs(r3(start+600:start+800)-uc(start+600:start+800))); 
% IAE_llm2_4=sum(abs(r3(start+800:1630)-uc(start+800:1630))); 
% 
% resid_1=r3(start:start+300)-uc(start:start+300);         count=find(resid_1>0);   IAE2_llm2_1=sum(resid_1(count));
% resid_2=r3(start+300:start+600)-uc(start+300:start+600); count2=find(resid_2>0);  IAE2_llm2_2=sum(resid_2(count2));
% resid_3=r3(start+600:start+800)-uc(start+600:start+800); count3=find(resid_3>0);  IAE2_llm2_3=sum(resid_3(count3));
% resid_4=r3(start+800:1630)-uc(start+800:1630);           count4=find(resid_4>0);  IAE2_llm2_4=sum(resid_4(count4));
% 
% load_llm2_1=sum(q_ramp(start:start+300))/6;
% load_llm2_2=sum(q_ramp(start+300:start+600))/6;
% load_llm2_3=sum(q_ramp(start+600:start+800))/6;
% load_llm2_4=sum(q_ramp(start+800:1630))/6;
% 
% indexs =find( (r3(start:start+300)<=uc(start+300)+0.01)     & (r3(start:start+300)>uc(start+300)-0.01) );     count=indexs-lag(indexs,1,0);     I=find(count==1);   Set_time_llm2_1=indexs(I(1)); 
% indexs2=find( (r3(start+300:start+600)<=uc(start+600)+0.01) & (r3(start+300:start+600)>uc(start+600)-0.01)); count2=indexs2-lag(indexs2,1,0);   I2=find(count2==1); Set_time_llm2_2=indexs2(I2(1));  
% indexs3=find( (r3(start+600:start+800)<=uc(start+800)+0.01) & (r3(start+600:start+800)>uc(start+800)-0.01)); count3=indexs3-lag(indexs3,1,0);   I3=find(count3==1); Set_time_llm2_3=indexs3(I3(1));
% indexs4=find( (r3(start+800:1630)<=uc(1630)+0.01)           & (r3(start+800:1630)<=uc(1630)-0.01));          count4=indexs4-lag(indexs4,1,0);   I4=find(count4==1); Set_time_llm2_4=indexs4(I4(1));
% 


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
          
            e(n) = uc(n) - r3(n);
            
            q_ramp(n) = q_ramp(n-1) + 15*e(n);
            
           if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)
        
         q_ramp(n)=min((p3(4)/r3(n))*10,20); if q_ramp(n)>20; q_ramp(n)=20; end   
    end

    waitbar(n/length(q1),h);
    
end

close(h)

% figure(3); subplot(211);h1=plot(r3,'m'); subplot(212);plot(q_ramp,'m');

IAE_alinea_1=sum(abs(r3(start:start+300)-uc(start:start+300))); 
IAE_alinea_2=sum(abs(r3(start+300:start+600)-uc(start+300:start+600))); 
IAE_alinea_3=sum(abs(r3(start+600:start+800)-uc(start+600:start+800))); 
IAE_alinea_4=sum(abs(r3(start+800:1630)-uc(start+800:1630))); 

resid_1=r3(start:start+300)-uc(start:start+300);         count=find(resid_1>0);   IAE2_alinea_1=sum(resid_1(count));
resid_2=r3(start+300:start+600)-uc(start+300:start+600); count2=find(resid_2>0);  IAE2_alinea_2=sum(resid_2(count2));
resid_3=r3(start+600:start+800)-uc(start+600:start+800); count3=find(resid_3>0);  IAE2_alinea_3=sum(resid_3(count3));
resid_4=r3(start+800:1630)-uc(start+800:1630);           count4=find(resid_4>0);  IAE2_alinea_4=sum(resid_4(count4));

load_alinea_1=sum(q_ramp(start:start+300))/6;
load_alinea_2=sum(q_ramp(start+300:start+600))/6;
load_alinea_3=sum(q_ramp(start+600:start+800))/6;
load_alinea_4=sum(q_ramp(start+800:1630))/6;

% indexs =find( (r3(start:start+300)<=uc(start+300)+0.01)     & (r3(start:start+300)>uc(start+300)-0.01) );    count=indexs-lag(indexs,1,0);     I=find(count==1);   Set_time_alinea_1=indexs(I(1));
% indexs2=find( (r3(start+300:start+600)<=uc(start+600)+0.01) & (r3(start+300:start+600)>uc(start+600)-0.01)); count2=indexs2-lag(indexs2,1,0);  I2=find(count2==1); Set_time_alinea_2=indexs2(I2(1));
% indexs3=find( (r3(start+600:start+800)<=uc(start+800)+0.01) & (r3(start+600:start+800)>uc(start+800)-0.01)); count3=indexs3-lag(indexs3,1,0);  I3=find(count3==1); Set_time_alinea_3=indexs3(I3(1));
% indexs4=find( (r3(start+800:1630)<=uc(1630)+0.01)           & (r3(start+800:1630)<=uc(1630)-0.01));          count4=indexs4-lag(indexs4,1,0);  I4=find(count4==1); Set_time_alinea_4=indexs4(I4(1));

load sdp_results

%(1) PARAMETRIC
par_A=zeros(size(q1)); par_B=par_A; par_C=par_A;numMFs=[2:14]; e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); 
% parametric_sdp
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
            
            F = [0 0 0 0;
                 par_B(n-2)  par_A(n-2)  par_C(n-2) 0;
                0 0 0 0;
                -par_B(n-2) -par_A(n-2) -par_C(n-2) 1];
            
            g =[0 pars(1,4) 0 -pars(1,4)];
            
            Q = eye(size(F)); Q(end,end)=250; r=10;  k=dlqri(F,g',Q,r); 
            
            k_1(i)=k(1); k_2(i) = k(2);  k_3(i) = k(3);  k_4(i) = -k(4); 
            
            e(n) = uc(n) - r3(n);
            
            q_ramp(n)=q_ramp(n-1) + k_4(i)*e(n) - k_2(i)*(r3(n)-r3(n-1)) - k_1(i)*(r2(n)-r2(n-1)) - k_3(i)*(r4(n)-r4(n-1));  
            
           if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
            
        else   
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)
        
         q_ramp(n)=min((p3(4)/r3(n))*10,20);   if q_ramp(n)>20; q_ramp(n)=20; end   
    end

    waitbar(n/length(q1),h);
    
end

close(h)

figure(3); 
subplot(311);h1=plot(r3(1:2100),'r'); grid
set(gca,'xtick',[1 600 1200 1800 2100], 'xticklabel',['06:00';'07:00';'08:00';'09:00'; '09:30']) 
axis([1 2100 0 150]);

subplot(312);h1=plot(r3(1:2100),'r'); axis([600 1800 40 70]); grid
set(gca,'xtick',[1 600 1200 1800 2100], 'xticklabel',['06:00';'07:00';'08:00';'09:00'; '09:30']) 

subplot(313);plot(q_ramp(1:2100),'r');  grid
set(gca,'xtick',[1 600 1200 1800 2100], 'xticklabel',['06:00';'07:00';'08:00';'09:00'; '09:30']) 
axis([1 2100 0 30]);

% uc=[zeros(600,1); (p3(4)-5)*ones(300,1); (p3(4)+5)*ones(300,1); (p3(4)+10)*ones(200,1); (p3(4)+2)*ones(800,1)];
% set(gca,'xtick',[1 50 100 150 200 250], 'xticklabel', ['07:40'; '08:10'; '08:40'; '09:10'; '09:40'; '10:10';]);grid
break
start=600; finish=2100;


IAE_sdp_1=sum(abs(r3(start:start+300)-uc(start:start+300))); 
IAE_sdp_2=sum(abs(r3(start+300:start+600)-uc(start+300:start+600))); 
IAE_sdp_3=sum(abs(r3(start+600:start+800)-uc(start+600:start+800))); 
IAE_sdp_4=sum(abs(r3(start+800:1630)-uc(start+800:1630))); 

resid_1=r3(start:start+300)-uc(start:start+300);         count=find(resid_1>0);   IAE2_sdp_1=sum(resid_1(count));
resid_2=r3(start+300:start+600)-uc(start+300:start+600); count2=find(resid_2>0);  IAE2_sdp_2=sum(resid_2(count2));
resid_3=r3(start+600:start+800)-uc(start+600:start+800); count3=find(resid_3>0);  IAE2_sdp_3=sum(resid_3(count3));
resid_4=r3(start+800:1630)-uc(start+800:1630);           count4=find(resid_4>0);  IAE2_sdp_4=sum(resid_4(count4));

load_sdp_1=sum(q_ramp(start:start+300))/6;
load_sdp_2=sum(q_ramp(start+300:start+600))/6;
load_sdp_3=sum(q_ramp(start+600:start+800))/6;
load_sdp_4=sum(q_ramp(start+800:1630))/6;

% indexs =find( (r3(start:start+300)<=uc(start+300)+0.01)     & (r3(start:start+300)>uc(start+300)-0.01) );    count=indexs-lag(indexs,1,0);     I=find(count==1);  Set_time_sdp_1=indexs(I(1)); 
% indexs2=find( (r3(start+300:start+600)<=uc(start+600)+0.01) & (r3(start+300:start+600)>uc(start+600)-0.01)); count2=indexs2-lag(indexs2,1,0); I2=find(count2==1); Set_time_sdp_2=indexs2(I2(1));
% indexs3=find( (r3(start+600:start+800)<=uc(start+800)+0.01) & (r3(start+600:start+800)>uc(start+800)-0.01)); count3=indexs3-lag(indexs3,1,0); I3=find(count3==1); Set_time_sdp_3=indexs3(I3(1));
% indexs4=find( (r3(start+800:1630)<=uc(1630)+0.01)           & (r3(start+800:1630)<=uc(1630)-0.01));          count4=indexs4-lag(indexs4,1,0); I4=find(count4==1); Set_time_sdp_4=indexs4(I4(1));

figure(5); 
o1=occ2den(r1,4,2,3,0); subplot(321); plot(q1,o1,'g.'); o2=occ2den(r2,4,2,3,0);     subplot(322); plot(q2,o2,'g.');
o3=occ2den(r3,4,2,3,0); subplot(323); plot(q3,o3,'g.'); o4=occ2den(r4,4,2,3,0);     subplot(324); plot(q4,o4,'g.');
o5=occ2den(r5,4,2,3,0); subplot(325); plot(q5,o5,'g.'); o6=occ2den(r_ramp,4,2,3,0); subplot(326); plot(q_ramp,o6,'g.');

figure(7);
subplot(311); plot(x1_irw,pars(:,1));hold; plot(state_A,par_A,'ro');
subplot(312); plot(x2_irw,pars(:,2));hold; plot(state_B,par_B,'ro');
subplot(313); plot(x3_irw,pars(:,3));hold; plot(state_C,par_C,'ro');

figure(8);
subplot(221);h=plot(k_1(2:end-10),'k'); set(h,'LineWidth',2);
set(gca,'xtick',[1 60 120 180 240], 'xticklabel', ['07:00'; '07:30'; '08:00'; '08:30'; '09:00';]);grid
h=ylabel('K_s_-_1'); set(h,'FontSize',14); h=xlabel('Time (hr)'); set(h,'FontSize',14);
subplot(222);h=plot(k_2(2:end-10),'k'); set(h,'LineWidth',2); 
set(gca,'xtick',[1 60 120 180 240], 'xticklabel', ['07:00'; '07:30'; '08:00'; '08:30'; '09:00';]);grid
h=ylabel('K_s'); set(h,'FontSize',14);h=xlabel('Time (hr)'); set(h,'FontSize',14);
subplot(223);h=plot(k_3(2:end-10),'k'); set(h,'LineWidth',2); 
set(gca,'xtick',[1 60 120 180 240], 'xticklabel', ['07:00'; '07:30'; '08:00'; '08:30'; '09:00';]);grid
h=ylabel('K_s_+_1'); set(h,'FontSize',14);h=xlabel('Time (hr)'); set(h,'FontSize',14);
subplot(224);h=plot(k_4(2:end-10),'k'); set(h,'LineWidth',2); 
set(gca,'xtick',[1 60 120 180 240], 'xticklabel', ['07:00'; '07:30'; '08:00'; '08:30'; '09:00';]);grid
h=ylabel('K_I'); set(h,'FontSize',14);h=xlabel('Time (hr)'); set(h,'FontSize',14);

% 49 veh/km (from 07:40 until 08:30)
% 59 veh/km (from 08:31 until 09:40)
% 64 veh/km (from 09:41 until 09:50)
% 56 veh/km (from 09:51 until 10:28)


break
%(2) Non-PARAMETRIC
e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); 
par_A_h=zeros(size(q1)); par_B_h=par_A_h; par_C_h=par_A_h; 
x1s_irw=irwsm(x1_irw,1,0.1); x2s_irw=irwsm(x2_irw,1,0.1); x3s_irw=irwsm(x3_irw,1,0.1);

h=waitbar(0,'Please wait.... Control in Progress (non-parametric) [5/7]');

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
                       
            % read-off States 
            state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x1s_irw); state_A(n-2)=min(x1s_irw); elseif state_A(n-2)>max(x1s_irw); state_A(n-2)=max(x1s_irw);end 
            state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x2s_irw); state_B(n-2)=min(x2s_irw); elseif state_B(n-2)>max(x2s_irw); state_B(n-2)=max(x2s_irw);end
            state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x3s_irw); state_C(n-2)=min(x3s_irw); elseif state_C(n-2)>max(x3s_irw); state_C(n-2)=max(x3s_irw);end
            
            % estimated non-parametric values
            par_A_h(n-2) =  interp1(x1s_irw,pars(:,1), state_A(n-2));
            par_B_h(n-2) =  interp1(x2s_irw,pars(:,2), state_B(n-2));
            par_C_h(n-2) =  interp1(x3s_irw,pars(:,3), state_C(n-2));
            
            F = [1 0 0 0;
                par_B_h(n-2) par_A_h(n-2) par_C_h(n-2) 0;
                0 0 1 0;
                -par_B_h(n-2) -par_A_h(n-2) -par_C_h(n-2) 1];
            
            g =[0 pars(1,4) 0 -pars(1,4)];
            
            Q = eye(size(F)); Q(end,end)=300; r=1;  k=dlqri(F,g',Q,r); 
            
            k_1(i)=k(1,1); k_2(i) = k(1,2);  k_3(i) = k(1,3);  k_4(i) = -k(1,4); 
            
            e(n) = uc(n) - r3(n);
            
            q_ramp(n)=q_ramp(n-1) + k_4(i)*e(n) - k_2(i)*(r3(n)-r3(n-1)) - k_1(i)*(r2(n)-r2(n-1)) - k_3(i)*(r4(n)-r4(n-1));  
            
           if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
            
        else
            q_ramp(n)=q_ramp(n-1);
        end
    end  
    if  n>=horizon(time)
        
         q_ramp(n)=min((p3(4)/r3(n))*10,20);   if q_ramp(n)>20; q_ramp(n)=20; end   
    end

    waitbar(n/length(q1),h);
    
end

close(h)

figure(3); subplot(211);h1=plot(r3,'c'); subplot(212);plot(q_ramp,'c');

IAE_sdp2_1=sum(abs(r3(start:start+300)-uc(start:start+300))); 
IAE_sdp2_2=sum(abs(r3(start+300:start+600)-uc(start+300:start+600))); 
IAE_sdp2_3=sum(abs(r3(start+600:start+800)-uc(start+600:start+800))); 
IAE_sdp2_4=sum(abs(r3(start+800:1630)-uc(start+800:1630))); 

resid_1=r3(start:start+300)-uc(start:start+300);         count=find(resid_1>0);   IAE2_sdp2_1=sum(resid_1(count));
resid_2=r3(start+300:start+600)-uc(start+300:start+600); count2=find(resid_2>0);  IAE2_sdp2_2=sum(resid_2(count2));
resid_3=r3(start+600:start+800)-uc(start+600:start+800); count3=find(resid_3>0);  IAE2_sdp2_3=sum(resid_3(count3));
resid_4=r3(start+800:1630)-uc(start+800:1630);           count4=find(resid_4>0);  IAE2_sdp2_4=sum(resid_4(count4));

load_sdp2_1=sum(q_ramp(start:start+300))/6;
load_sdp2_2=sum(q_ramp(start+300:start+600))/6;
load_sdp2_3=sum(q_ramp(start+600:start+800))/6;
load_sdp2_4=sum(q_ramp(start+800:1630))/6;

% indexs =find( (r3(start:start+300)<=uc(start+300)+0.01)     & (r3(start:start+300)>uc(start+300)-0.01) );    count=indexs-lag(indexs,1,0);     I=find(count==1);  Set_time_sdp2_1=indexs(I(1));
% indexs2=find( (r3(start+300:start+600)<=uc(start+600)+0.01) & (r3(start+300:start+600)>uc(start+600)-0.01)); count2=indexs2-lag(indexs2,1,0); I2=find(count2==1); Set_time_sdp2_2=indexs2(I2(1));
% indexs3=find( (r3(start+600:start+800)<=uc(start+800)+0.01) & (r3(start+600:start+800)>uc(start+800)-0.01)); count3=indexs3-lag(indexs3,1,0); I3=find(count3==1); Set_time_sdp2_3=indexs3(I3(1));
% indexs4=find( (r3(start+800:1630)<=uc(1630)+0.01)           & (r3(start+800:1630)<=uc(1630)-0.01));          count4=indexs4-lag(indexs4,1,0); I4=find(count4==1); Set_time_sdp2_4=indexs4(I4(1));
% 


% %(2) STM-PIP
% e=zeros(length(q1),1); state_A=zeros(size(q1)); state_B=zeros(size(q1)); state_C=zeros(size(q1)); 
% % par_A_h=zeros(size(q1)); par_B_h=par_A_h; par_C_h=par_A_h; 
% % x1s_irw=irwsm(x1_irw,1,0.1); x2s_irw=irwsm(x2_irw,1,0.1); x3s_irw=irwsm(x3_irw,1,0.1);
% 
% h=waitbar(0,'Please wait.... Control in Progress (non-parametric) [7/7]');
% 
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
%             i=i+1;
%             e(n) = uc(n) - r3(n);
%                        
% %             % read-off States 
% %             state_A(n-2) = max(r3(n-2),r4(n-2)); if state_A(n-2) <= min(x1s_irw); state_A(n-2)=min(x1s_irw); elseif state_A(n-2)>max(x1s_irw); state_A(n-2)=max(x1s_irw);end 
% %             state_B(n-2) = max(r2(n-2),r3(n-2)); if state_B(n-2) <= min(x2s_irw); state_B(n-2)=min(x2s_irw); elseif state_B(n-2)>max(x2s_irw); state_B(n-2)=max(x2s_irw);end
% %             state_C(n-2) = max(r4(n-2),r5(n-2)); if state_C(n-2) <= min(x3s_irw); state_C(n-2)=min(x3s_irw); elseif state_C(n-2)>max(x3s_irw); state_C(n-2)=max(x3s_irw);end
% %             
% %             % estimated non-parametric values
% %             par_A_h(n-2) =  interp1(x1s_irw,pars(:,1), state_A(n-2));
% %             par_B_h(n-2) =  interp1(x2s_irw,pars(:,2), state_B(n-2));
% %             par_C_h(n-2) =  interp1(x3s_irw,pars(:,3), state_C(n-2));
% %             
% %             F = [0 0 0 0;
% %                 par_B_h(n-2) par_A_h(n-2) par_C_h(n-2) 0;
% %                 0 0 0 0;
% %                 -par_B_h(n-2) -par_A_h(n-2) -par_C_h(n-2) 1];
% %             
% %             g =[0 pars(1,4) 0 -pars(1,4)];
% %             
% %             Q = eye(size(F)); Q(end,end)=100; r=1;  k=dlqri(F,g',Q,r); 
% %             
% %             k_1(n)=k(1,1); k_2(n) = k(1,2);  k_3(n) = k(1,3);  k_4(n) = -k(1,4); 
% %             
% %             e(n) = uc(n) - r3(n);  
% %             
% %             q_ramp(n)=q_ramp(n-1) + k_4(n)*e(n) - k_2(n)*(r3(n)-r3(n-1)) - k_1(n)*(r2(n)-r2(n-1)) - k_3(n)*(r4(n)-r4(n-1));  
%             
%            if q_ramp(n)<0; q_ramp(n)=0; elseif q_ramp(n)>20; q_ramp(n)=20;  end  
%             
%         else   
%             q_ramp(n)=q_ramp(n-1);
%         end
%     end  
%     if  n>=horizon(time)
%         
%          q_ramp(n)=min((p3(4)/r3(n))*10,20);   if q_ramp(n)>20; q_ramp(n)=20; end   
%     end
% 
%     waitbar(n/length(q1),h);
%     
% end
% 
% close(h)
% 
% figure(3); subplot(211);h1=plot(r3,'c'); subplot(212);plot(q_ramp,'c');
% 
% IAE_stm_1=sum(abs(r3(start:start+300)-uc(start:start+300))); 
% IAE_stm_2=sum(abs(r3(start+300:start+600)-uc(start+300:start+600))); 
% IAE_stm_3=sum(abs(r3(start+600:start+800)-uc(start+600:start+800))); 
% IAE_stm_4=sum(abs(r3(start+800:1630)-uc(start+800:1630))); 
% 
% resid_1=r3(start:start+300)-uc(start:start+300);         count=find(resid_1>0);   IAE2_stm_1=sum(resid_1(count));
% resid_2=r3(start+300:start+600)-uc(start+300:start+600); count2=find(resid_2>0);  IAE2_stm_2=sum(resid_2(count2));
% resid_3=r3(start+600:start+800)-uc(start+600:start+800); count3=find(resid_3>0);  IAE2_stm_3=sum(resid_3(count3));
% resid_4=r3(start+800:1630)-uc(start+800:1630);           count4=find(resid_4>0);  IAE2_stm_4=sum(resid_4(count4));
% 
% load_stm_1=sum(q_ramp(start:start+300))/6;
% load_stm_2=sum(q_ramp(start+300:start+600))/6;
% load_stm_3=sum(q_ramp(start+600:start+800))/6;
% load_stm_4=sum(q_ramp(start+800:1630))/6;
% 
% indexs =find( (r3(start:start+300)<=uc(start+300)+0.01)     & (r3(start:start+300)>uc(start+300)-0.01) );     Set_time_stm_1=indexs(2);
% indexs2=find( (r3(start+300:start+600)<=uc(start+600)+0.01) & (r3(start+300:start+600)>uc(start+600)-0.01));  Set_time_stm_2=indexs2(2);
% indexs3=find( (r3(start+600:start+800)<=uc(start+800)+0.01) & (r3(start+600:start+800)>uc(start+800)-0.01));  Set_time_stm_3=indexs3(2);
% indexs4=find( (r3(start+800:1630)<=uc(1630)+0.01)           & (r3(start+800:1630)<=uc(1630)-0.01));           Set_time_stm_4=indexs4(2);
% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IAE_llm  = [IAE_llm_1 IAE_llm_2 IAE_llm_3 IAE_llm_4]';
IAE2_llm = [IAE2_llm_1 IAE2_llm_2 IAE2_llm_3 IAE2_llm_4]';
Load_llm = [load_llm_1 load_llm_2 load_llm_3 load_llm_4]';
Set_time_llm=[Set_time_llm_1 Set_time_llm_2 Set_time_llm_3 Set_time_llm_4]';

IAE_llm2  = [IAE_llm2_1 IAE_llm2_2 IAE_llm2_3 IAE_llm2_4]';
IAE2_llm2 = [IAE2_llm2_1 IAE2_llm2_2 IAE2_llm2_3 IAE2_llm2_4]';
Load_llm2 = [load_llm2_1 load_llm2_2 load_llm2_3 load_llm2_4]';
Set_time_llm2=[Set_time_llm2_1 Set_time_llm2_2 Set_time_llm2_3 Set_time_llm2_4]';

IAE_alinea  = [IAE_alinea_1 IAE_alinea_2 IAE_alinea_3 IAE_alinea_4]';
IAE2_alinea = [IAE2_alinea_1 IAE2_alinea_2 IAE2_alinea_3 IAE2_alinea_4]';
Load_alinea = [load_alinea_1 load_alinea_2 load_alinea_3 load_alinea_4]';
Set_time_alinea=[Set_time_alinea_1 Set_time_alinea_2 Set_time_alinea_3 Set_time_alinea_4]';

IAE_sdp  = [IAE_sdp_1 IAE_sdp_2 IAE_sdp_3 IAE_sdp_4]';
IAE2_sdp = [IAE2_sdp_1 IAE2_sdp_2 IAE2_sdp_3 IAE2_sdp_4]';
Load_sdp = [load_sdp_1 load_sdp_2 load_sdp_3 load_sdp_4]';
Set_time_sdp=[Set_time_sdp_1 Set_time_sdp_2 Set_time_sdp_3 Set_time_sdp_4]';

IAE_sdp2  = [IAE_sdp2_1 IAE_sdp2_2 IAE_sdp2_3 IAE_sdp2_4]';
IAE2_sdp2 = [IAE2_sdp2_1 IAE2_sdp2_2 IAE2_sdp2_3 IAE2_sdp2_4]';
Load_sdp2 = [load_sdp2_1 load_sdp2_2 load_sdp2_3 load_sdp2_4]';
Set_time_sdp2=[Set_time_sdp2_1 Set_time_sdp2_2 Set_time_sdp2_3 Set_time_sdp2_4]';

% IAE_stm  = [IAE_stm_1 IAE_stm_2 IAE_stm_3 IAE_stm_4]';
% IAE2_stm = [IAE2_stm_1 IAE2_stm_2 IAE2_stm_3 IAE2_stm_4]';
% Load_stm = [load_stm_1 load_stm_2 load_stm_3 load_stm_4]';
% Set_time_stm=[Set_time_stm_1 Set_time_stm_2 Set_time_stm_3 Set_time_stm_4]';

figure(1);
subplot(221);bar([IAE_llm IAE_alinea IAE_sdp IAE_sdp2 ]); h=ylabel('IAE'); set(h,'FontSize',14);
subplot(222);bar([IAE2_llm IAE2_alinea IAE2_sdp IAE2_sdp2 ]);h=ylabel('IAE +'); set(h,'FontSize',14);
subplot(223);bar([Load_llm Load_alinea Load_sdp  Load_sdp2]); h=ylabel('Control Effort'); set(h,'FontSize',14);
subplot(224);bar([Set_time_llm Set_time_alinea Set_time_sdp Set_time_sdp2]);
h=ylabel('Settling Time 1%'); set(h,'FontSize',14); legend('LLM-PIP','LLM-PIP2','ALINEA','SDP-PIP','SDP-PIP2')