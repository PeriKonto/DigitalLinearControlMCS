%function [f, g, k]=gains(Aq, Bq, v)
 function [f, g, k]=gains(a, b, v)
 
 % finds F, G and k gains from the control gain vector v
%
% b/a model (trunciated form)
% v from pip1.m or pipopt.m

% James Taylor
% Modified 28.9.94
% Original 21.9.94

lv=length(v);
la=length(a);

f=v(1:la)';
g=[1 v(la+1:lv-1)'];

k=v(lv);

%end
