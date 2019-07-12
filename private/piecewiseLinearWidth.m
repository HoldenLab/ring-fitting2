function [w] = piecewiseLinearWidth(a,t)
%piecewiseLinearWidth 
% fit flat line and then linear constriction
% t0 = a(1);
% t1 = a(2);
% d = a(3);

t0 = a(1);
t1 = a(2);
d = a(3);
w=zeros(size(t));

%split into 2 portions
%t<t0
w(t<t0)=d;
%t>=t0
w(t>=t0) = d*(1 - (t(t>=t0)-t0)/(t1-t0));

end

