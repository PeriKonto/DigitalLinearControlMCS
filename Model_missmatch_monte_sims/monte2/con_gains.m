function [f,g,k]=con_gains(A2,B2,at,bt)
%at =[1 -0.2146];
%bt = [0 14.5442];
%at
%bt
a = at(2:length(at));
b = bt(2:length(bt));


%optimal control procedure
Aq= [-A2 0; A2 1;];
Bq= [B2 -B2]';
%A2
%B2
%Aq = [0.2146 0;-0.2146 1;];                   
%Bq = [14.5442 -14.5442]'; %B is g in state-space form

D = [0 1];
R = 0.01;
h = [1  0];
Q = [0.01 0;0 1];
[v,S,E] = dlqr(Aq,Bq,Q,R);

%find the system model
F_closed = (Aq-Bq*v);
d_closed = zeros(size(h,1),size(D,2));
sys = ss(F_closed,D',h,0,-1);
P = pole(sys);
%figure(2)
%pzmap(sys);
%v=v';
%v = -v;
[f, g, k]=gains(a,b,v);