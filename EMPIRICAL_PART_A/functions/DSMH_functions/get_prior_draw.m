function [draw,priorval] = get_prior_draw(c,m,s,nu,signA)

draw=zeros(1,c);

for jx=1:c
    if signA(jx,1) == 0
       draw(1,jx)=draw_student(m(jx,1),s(jx,1),nu(jx,1));
    else
       draw(1,jx)=draw_truncated_student(m(jx,1),s(jx,1),nu(jx,1),signA(jx,1));
    end
end

priorval=exp(logP_factor(draw',m,s,nu,signA));
    