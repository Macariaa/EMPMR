function [errR, errT] = err_cal(grtR, grtT, pM, scannum, s)

errR = 0;
errT = 0;
for i=1:scannum
    p(i).M(1:4,1:3) = [grtR{1,i};0 0 0];
    p(i).M(1:4,4) = [grtT{1,i}*s;1];
    p(i).M= inv(p(1).M)* p(i).M;
    R1= pM(i).M(1:3,1:3);
    R2= p(i).M(1:3,1:3);
    R12= trace(R1'*R2);
    errR = errR + real(acos(0.5*(R12-1)));
    errT = errT + norm(pM(i).M(1:3,4)-p(i).M(1:3,4),2);
end
errR = errR/scannum;
errT = errT/scannum;

function q= R2vec(R)
t=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2;
q=[t (R(3,2)-R(2,3))/(4*t) (R(1,3)-R(3,1))/(4*t) (R(2,1)-R(1,2))/(4*t)];
