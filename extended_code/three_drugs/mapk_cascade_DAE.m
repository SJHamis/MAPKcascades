function yp = mapk_cascade_DAE(y, BRAF_in, ATP_in, DBF_in, TMT_in, SCH_in)

BRAF_tot=BRAF_in;
ATP_tot=ATP_in;
MEK_tot=1.2;
ERK_tot=1.2;
O_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
phosph3_tot=0.12;
DBF_tot=DBF_in;
TMT_tot=TMT_in;
SCH_tot=SCH_in;
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
a11=a_HF;
a12=a_HF;
d1=6.23; %0.8268; % ATP binding to BRAF (or BRAF complex).
d2=0.02385; % MEK binding to BRAF (or BRAF complex).
d5=0.3468; % ATP binding to MEK (or MEK complex).
d6=0.01184; % ERK binding to MEK (or MEK complex).
d4=0.00005936; % DBF binding to BRAF (or BRAF complex).
d8=0.0012296; % TMT binding to MEK (or MEK complex).
d3=0.0159; % BRAF phosphasate.
d7=0.0159; % MEK phosphasate.
d9=14.64; % ATP binding to ERK (or ERK complex).
d10=0.7738; % S binding to ERK (or ERK complex).
d11=0.0159; % ERK phosphatase.
d12=0.002967; % SCH binding to ERK (or ERK complex).
k12=0.66;
k56=0.0242;
k3=0.0159;
k7=0.0159;
k910=0.2;
k11=0.0159;

yp=[y(1)+y(3)+y(5)+y(7)+y(8)+y(10)+y(16)+y(17)+y(18)-BRAF_tot
y(2)+y(3)+y(8)+y(10)+y(19)+y(24)+y(25)+y(32)+y(35)+y(36)+y(37)+y(42)+y(43)-ATP_tot
+y(3)*(-d1)+a1*(y(1)*y(2))+y(3)*y(4)*(-a2)+d2*y(8)+y(3)*y(6)*(-a2)+d2*y(10)
y(4)+y(5)+y(6)+y(7)+y(8)+y(10)+y(11)+y(13)+y(14)+y(17)+y(18)+y(19)+y(21)+y(23)+y(24)+y(25)+y(31)+y(32)+y(33)+y(34)+y(35)+y(36)-MEK_tot
+y(5)*(-d2)+a2*(y(1)*y(4))+y(5)*y(2)*(-a1)+d1*y(8)+y(5)*y(15)*(-a4)+d4*y(17)
+y(1)*y(6)*(-a2)+d2*y(7)+y(8)*k12+y(3)*y(6)*(-a2)+d2*y(10)+y(8)*k12+y(6)*y(12)*(-a3)+d3*y(13)+y(14)*k3+y(16)*y(6)*(-a2)+d2*y(18)
+y(7)*(-d2)+a2*(y(1)*y(6))+y(7)*y(2)*(-a1)+d1*y(10)+y(7)*y(15)*(-a4)+d4*y(18)
+y(8)*(-d2-k12)+a2*(y(3)*y(4))+y(8)*(-d1-k12)+a1*(y(5)*y(2))
+y(8)*k12+y(10)*k12+y(8)*k12+y(10)*k12+y(24)*k56+y(25)*k56+y(24)*k56+y(25)*k56+y(42)*k910+y(43)*k910+y(42)*k910+y(43)*k910
+y(10)*(-d2-k12)+a2*(y(3)*y(6))+y(10)*(-d1-k12)+a1*(y(7)*y(2))
+y(10)*k12+y(10)*k12+y(11)*y(12)*(-a3)+d3*y(14)+y(11)*y(2)*(-a5)+d5*y(19)+y(11)*y(20)*(-a6)+d6*y(21)+y(11)*y(22)*(-a6)+d6*y(23)+y(24)*k56+y(25)*k56+y(24)*k56+y(25)*k56+y(11)*y(30)*(-a8)+d8*y(31)
y(12)+y(13)+y(14)-phosph1_tot
+y(13)*(-d3-k3)+a3*(y(6)*y(12))
+y(14)*(-d3-k3)+a3*(y(11)*y(12))
y(15)+y(16)+y(17)+y(18)-DBF_tot
+y(16)*(-d4)+a4*(y(1)*y(15))+y(16)*y(4)*(-a2)+d2*y(17)+y(16)*y(6)*(-a2)+d2*y(18)
+y(17)*(-d4)+a4*(y(5)*y(15))+y(17)*(-d2)+a2*(y(16)*y(4))
+y(18)*(-d4)+a4*(y(7)*y(15))+y(18)*(-d2)+a2*(y(16)*y(6))
+y(19)*(-d5)+a5*(y(11)*y(2))+y(19)*y(20)*(-a6)+d6*y(24)+y(19)*y(22)*(-a6)+d6*y(25)+y(19)*y(30)*(-a8)+d8*y(32)
y(20)+y(21)+y(22)+y(23)+y(24)+y(25)+y(26)+y(28)+y(29)+y(33)+y(34)+y(35)+y(36)+y(37)+y(39)+y(41)+y(42)+y(43)+y(49)+y(50)+y(51)-ERK_tot
+y(21)*(-d6)+a6*(y(11)*y(20))+y(21)*y(2)*(-a5)+d5*y(24)+y(21)*y(30)*(-a8)+d8*y(33)
+y(11)*y(22)*(-a6)+d6*y(23)+y(24)*k56+y(19)*y(22)*(-a6)+d6*y(25)+y(24)*k56+y(22)*y(27)*(-a7)+d7*y(28)+y(29)*k7+y(31)*y(22)*(-a6)+d6*y(34)+y(32)*y(22)*(-a6)+d6*y(36)
+y(23)*(-d6)+a6*(y(11)*y(22))+y(23)*y(2)*(-a5)+d5*y(25)+y(23)*y(30)*(-a8)+d8*y(34)
+y(24)*(-d6-k56)+a6*(y(19)*y(20))+y(24)*(-d5-k56)+a5*(y(21)*y(2))+y(24)*y(30)*(-a8)+d8*y(35)
+y(25)*(-d6-k56)+a6*(y(19)*y(22))+y(25)*(-d5-k56)+a5*(y(23)*y(2))+y(25)*y(30)*(-a8)+d8*y(36)
+y(25)*k56+y(25)*k56+y(26)*y(27)*(-a7)+d7*y(29)+y(26)*y(2)*(-a9)+d9*y(37)+y(26)*y(38)*(-a10)+d10*y(39)+y(26)*y(40)*(-a10)+d10*y(41)+y(42)*k910+y(43)*k910+y(42)*k910+y(43)*k910+y(26)*y(48)*(-a12)+d12*y(49)
y(27)+y(28)+y(29)-phosph2_tot
+y(28)*(-d7-k7)+a7*(y(22)*y(27))
+y(29)*(-d7-k7)+a7*(y(26)*y(27))
y(30)+y(31)+y(32)+y(33)+y(34)+y(35)+y(36)-TMT_tot
+y(31)*(-d8)+a8*(y(11)*y(30))+y(31)*y(2)*(-a5)+d5*y(32)+y(31)*y(20)*(-a6)+d6*y(33)+y(31)*y(22)*(-a6)+d6*y(34)
+y(32)*(-d8)+a8*(y(19)*y(30))+y(32)*(-d5)+a5*(y(31)*y(2))+y(32)*y(20)*(-a6)+d6*y(35)+y(32)*y(22)*(-a6)+d6*y(36)
+y(33)*(-d8)+a8*(y(21)*y(30))+y(33)*(-d6)+a6*(y(31)*y(20))+y(33)*y(2)*(-a5)+d5*y(35)
+y(34)*(-d8)+a8*(y(23)*y(30))+y(34)*(-d6)+a6*(y(31)*y(22))+y(34)*y(2)*(-a5)+d5*y(36)
+y(35)*(-d8)+a8*(y(24)*y(30))+y(35)*(-d5)+a5*(y(33)*y(2))+y(35)*(-d6)+a6*(y(32)*y(20))
+y(36)*(-d8)+a8*(y(25)*y(30))+y(36)*(-d5)+a5*(y(34)*y(2))+y(36)*(-d6)+a6*(y(32)*y(22))
+y(37)*(-d9)+a9*(y(26)*y(2))+y(37)*y(38)*(-a10)+d10*y(42)+y(37)*y(40)*(-a10)+d10*y(43)
y(38)+y(39)+y(40)+y(41)+y(42)+y(43)+y(44)+y(46)+y(47)+y(50)+y(51)-O_tot
+y(39)*(-d10)+a10*(y(26)*y(38))+y(39)*y(2)*(-a9)+d9*y(42)+y(39)*y(48)*(-a12)+d12*y(50)
+y(26)*y(40)*(-a10)+d10*y(41)+y(42)*k910+y(37)*y(40)*(-a10)+d10*y(43)+y(42)*k910+y(40)*y(45)*(-a11)+d11*y(46)+y(47)*k11+y(49)*y(40)*(-a10)+d10*y(51)
+y(41)*(-d10)+a10*(y(26)*y(40))+y(41)*y(2)*(-a9)+d9*y(43)+y(41)*y(48)*(-a12)+d12*y(51)
+y(42)*(-d10-k910)+a10*(y(37)*y(38))+y(42)*(-d9-k910)+a9*(y(39)*y(2))
+y(43)*(-d10-k910)+a10*(y(37)*y(40))+y(43)*(-d9-k910)+a9*(y(41)*y(2))
+y(43)*k910+y(43)*k910+y(44)*y(45)*(-a11)+d11*y(47)
y(45)+y(46)+y(47)-phosph3_tot
+y(46)*(-d11-k11)+a11*(y(40)*y(45))
+y(47)*(-d11-k11)+a11*(y(44)*y(45))
y(48)+y(49)+y(50)+y(51)-SCH_tot
+y(49)*(-d12)+a12*(y(26)*y(48))+y(49)*y(38)*(-a10)+d10*y(50)+y(49)*y(40)*(-a10)+d10*y(51)
+y(50)*(-d12)+a12*(y(39)*y(48))+y(50)*(-d10)+a10*(y(49)*y(38))
+y(51)*(-d12)+a12*(y(41)*y(48))+y(51)*(-d10)+a10*(y(49)*y(40))];

end