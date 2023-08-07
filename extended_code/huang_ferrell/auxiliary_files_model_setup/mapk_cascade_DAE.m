function yp = mapk_cascade_DAE(y, KKK_in, E1_in, E2_in)

KKK_tot=KKK_in;
E1_tot=E1_in;
E2_tot=E2_in;
KK_tot=1.2;
K_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
a_HF=0.106; %All a-values are the same.
a1=a_HF;
a2=a_HF;
a3=a_HF;
a4=a_HF;
a5=a_HF;
a6=a_HF;
a7=a_HF;
a8=a_HF;
a9=a_HF;
a10=a_HF;
d1=6.23; %0.8268; % E1 and E2 binding to KKK (or KKK complex).
d2=d1;
d3=0.02385; % KK binding to KKK (or KKK complex).
d5=d3;
d7=0.01184; % K binding to KK (or KK complex).
d9=d7;
d4=0.0159; % KKK phosphasate.
d6=d4;
d8=0.0159; % KK phosphasate.
d10=d8;
k1=0.66;
k2=k1;
k3=k1;
k5=k1;
k7=0.0242;
k9=k7;
k4=0.0159;
k6=k4;
k8=0.0159;
k10=k8;

yp=[+y(1)*y(2)*(-a1)+d1*y(3)+y(6)*k2
+y(1)*y(2)*(-a1)+d1*y(3)+y(3)*k1
+y(3)*(-d1-k1)+a1*(y(1)*y(2))
+y(3)*k1+y(4)*y(5)*(-a2)+d2*y(6)+y(7)*y(4)*(-a3)+d3*y(8)+y(8)*k3+y(9)*y(4)*(-a5)+d5*y(10)+y(10)*k5
+y(4)*y(5)*(-a2)+d2*y(6)+y(6)*k2
+y(6)*(-d2-k2)+a2*(y(4)*y(5))
+y(7)*y(4)*(-a3)+d3*y(8)+y(13)*k4
+y(8)*(-d3-k3)+a3*(y(7)*y(4))
+y(8)*k3+y(9)*y(4)*(-a5)+d5*y(10)+y(9)*y(12)*(-a4)+d4*y(13)+y(14)*k6
+y(10)*(-d5-k5)+a5*(y(9)*y(4))
+y(10)*k5+y(11)*y(12)*(-a6)+d6*y(14)+y(11)*y(15)*(-a7)+d7*y(16)+y(16)*k7+y(17)*y(11)*(-a9)+d9*y(18)+y(18)*k9
+y(9)*y(12)*(-a4)+d4*y(13)+y(13)*k4+y(11)*y(12)*(-a6)+d6*y(14)+y(14)*k6
+y(13)*(-d4-k4)+a4*(y(9)*y(12))
+y(14)*(-d6-k6)+a6*(y(11)*y(12))
+y(11)*y(15)*(-a7)+d7*y(16)+y(21)*k8
+y(16)*(-d7-k7)+a7*(y(11)*y(15))
+y(16)*k7+y(17)*y(11)*(-a9)+d9*y(18)+y(17)*y(20)*(-a8)+d8*y(21)+y(22)*k10
+y(18)*(-d9-k9)+a9*(y(17)*y(11))
+y(18)*k9+y(19)*y(20)*(-a10)+d10*y(22)
+y(17)*y(20)*(-a8)+d8*y(21)+y(21)*k8+y(19)*y(20)*(-a10)+d10*y(22)+y(22)*k10
+y(21)*(-d8-k8)+a8*(y(17)*y(20))
+y(22)*(-d10-k10)+a10*(y(19)*y(20))];

end