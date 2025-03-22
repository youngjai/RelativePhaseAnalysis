%
% Code for the solution of the
% differential equation:
%
%   dy/dt = f(t,y)
%


function [t,y] = rk2(f, init, T, n)

%  inputs

y(:,1)=init;
n =n -1;
 h = (max(T)-min(T))/n;

t = min(T):h:max(T);
%
%  rk2 loop
%
for j = 1:n
   k1=h * f(t(j),y(:,j));
   k2=h * f(t(j)+h/2,y(:,j)+k1/2);
   y(:,j+1)=y(:,j)+k2;
end

y=y';
t=t';