function y=NewSubs(f,x,val)
%
%
y=subs(f,x,val);

if ~prod(double(size(y)==size(val)))
    y=y*ones(size(val));
end
end