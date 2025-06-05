function m = CiMultiple(alpha_s,alpha0,alpha1)

a=(alpha0+alpha1)/2;
a0=a^5;
a1=a^.2;

if alpha_s <= a0
    m = 0.2;
elseif alpha_s >= a1
    m = 5;
else 
    m = log(a)/log(alpha_s);
end