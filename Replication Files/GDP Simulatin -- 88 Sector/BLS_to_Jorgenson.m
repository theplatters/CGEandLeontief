%% BLS IO DATA TO JORGENSON SECTORS 2006-2008

function [grossalt, vaddalt, dweightalt]=BLS_to_Jorgenson(grossy, vadd);


%% 1) READ BLS IO DATA
%Loading and setting up IO matrices for 1993-2008 BLS IO Matrices
% LOAD RAW DATA
load NOMINAL_USE_1993.DAT;
io(:,:,1)=NOMINAL_USE_1993;
load NOMINAL_USE_1994.DAT;
io(:,:,2)=NOMINAL_USE_1994;
load NOMINAL_USE_1995.DAT;
io(:,:,3)=NOMINAL_USE_1995;
load NOMINAL_USE_1996.DAT;
io(:,:,4)=NOMINAL_USE_1996;
load NOMINAL_USE_1997.DAT;
io(:,:,5)=NOMINAL_USE_1997;
load NOMINAL_USE_1998.DAT;
io(:,:,6)=NOMINAL_USE_1998;
load NOMINAL_USE_1999.DAT;
io(:,:,7)=NOMINAL_USE_1999;
load NOMINAL_USE_2000.DAT;
io(:,:,8)=NOMINAL_USE_2000;
load NOMINAL_USE_2001.DAT;
io(:,:,9)=NOMINAL_USE_2001;
load NOMINAL_USE_2002.DAT;
io(:,:,10)=NOMINAL_USE_2002;
load NOMINAL_USE_2003.DAT;
io(:,:,11)=NOMINAL_USE_2003;
load NOMINAL_USE_2004.DAT;
io(:,:,12)=NOMINAL_USE_2004;
load NOMINAL_USE_2005.DAT;
io(:,:,13)=NOMINAL_USE_2005;
load NOMINAL_USE_2006.DAT;
io(:,:,14)=NOMINAL_USE_2006;
load NOMINAL_USE_2007.DAT;
io(:,:,15)=NOMINAL_USE_2007;
load NOMINAL_USE_2008.DAT;
io(:,:,16)=NOMINAL_USE_2008;

clear NOMINAL_USE_1993 NOMINAL_USE_1994 NOMINAL_USE_1995 NOMINAL_USE_1996 NOMINAL_USE_1997 NOMINAL_USE_1998 NOMINAL_USE_1999 NOMINAL_USE_2000;
clear NOMINAL_USE_2001 NOMINAL_USE_2002 NOMINAL_USE_2003 NOMINAL_USE_2004 NOMINAL_USE_2005 NOMINAL_USE_2006 NOMINAL_USE_2007 NOMINAL_USE_2008;

%% 2) AGGREGATE BLS DATA TO JORGENSON 88 INDUSTRIES
% Aggregate BLS 202 industries data to Jorgenson 88 industries (for every
% year in the BLS sector data)
% Correspondence weights supplied by Dale Jorgenson and Associates

blstodj=zeros(203,203,16);
% Aggregating Rows of BLS IO data
for k=1:16
blstodj(:,1,k)=sum(io(:,1:2,k)')';
blstodj(:,2,k)=(io(:,6,k));
blstodj(:,3,k)=(io(:,3,k))+(io(:,5,k));
blstodj(:,4,k)=0.16*sum(io(:,8:10,k)')';
blstodj(:,5,k)=0.475*sum(io(:,8:10,k)')';
blstodj(:,6,k)=0.336*sum(io(:,8:10,k)')'+0.016*(io(:,11,k));
blstodj(:,7,k)=(io(:,7,k))+0.971*(io(:,11,k));
%blstodj(:,8,k)= uranium ore mining; no data in original Jorgenson database
blstodj(:,9,k)=0.235*0.915*(io(:,15,k))+0.109*(io(:,145,k))+0.235*0.054*(io(:,15,k))+0.134*0.83*(io(:,15,k))+0.134*0.054*(io(:,15,k))+0.632*0.984*(io(:,15,k));
blstodj(:,10,k)=sum(io(:,37:39,k)')'+(io(:,4,k));
blstodj(:,11,k)=sum(io(:,95:97,k)')';
blstodj(:,12,k)=sum(io(:,53:56,k)')';
blstodj(:,13,k)=sum(io(:,57:61,k)')';
blstodj(:,14,k)=sum(io(:,62:70,k)')';
blstodj(:,15,k)=sum(io(:,71:72,k)')'+sum(io(:,74:77,k)')'+(0.81-0.159)*io(:,73,k);
blstodj(:,16,k)=io(:,78,k);
blstodj(:,17,k)=0.018*(io(:,63,k))+0.975*(io(:,84,k))+0.399*(io(:,87,k))+0.023*(io(:,90,k));
blstodj(:,18,k)=io(:,80,k)+0.078*0.211*(io(:,113,k));
blstodj(:,19,k)=sum(io(:,85:87,k)')'+0.99*(io(:,84,k));
blstodj(:,20,k)=(io(:,79,k));
blstodj(:,21,k)=(io(:,81,k))+(io(:,83,k))+0.01*(io(:,136,k));
blstodj(:,22,k)=sum(io(:,88:90,k)')'+0.016*(io(:,136,k));
blstodj(:,23,k)=(io(:,91,k));
blstodj(:,24,k)=(io(:,93,k));
blstodj(:,25,k)=(io(:,92,k))+(io(:,94,k));
blstodj(:,26,k)=io(:,82,k)+0.19*(io(:,73,k));
blstodj(:,27,k)=(io(:,98,k));
blstodj(:,28,k)=0.314*(io(:,50,k))+0.159*(io(:,73,k))+0.01*(io(:,84,k));
blstodj(:,29,k)=0.899*(io(:,98,k));
blstodj(:,30,k)=sum(io(:,16:25,k)')'+0.031*(io(:,136,k));
blstodj(:,31,k)=(io(:,26,k));
blstodj(:,32,k)=sum(io(:,27:32,k)')';
blstodj(:,33,k)=0.982*(io(:,33,k))+0.944*(io(:,34,k));
blstodj(:,34,k)=0.056*(io(:,34,k))+0.982*(io(:,35,k))+0.89*(io(:,36,k));
blstodj(:,35,k)=sum(io(:,40:41,k)')'+0.013*(io(:,136,k));
blstodj(:,36,k)=(io(:,111,k))+0.078*0.126*(io(:,113,k))+(1-0.699-0.07-0.057-0.048)*(io(:,42,k));
blstodj(:,37,k)=(0.699+0.07+0.057+0.048)*(io(:,42,k));
blstodj(:,38,k)=sum(io(:,44:46,k)')'+sum(io(:,48:49,k)')'+0.458*(io(:,50,k))+0.014*(io(:,136,k));
blstodj(:,39,k)=io(:,47,k)+0.021*(io(:,133,k));
blstodj(:,40,k)=(io(:,43,k));
blstodj(:,41,k)=sum(io(:,51:52,k)')'+0.11*(io(:,36,k))+0.101*(io(:,99,k));
blstodj(:,42,k)=io(:,103,k)+0.102*0.11*(io(:,108,k));
blstodj(:,43,k)=io(:,106,k)+0.118*0.032*(io(:,108,k));
blstodj(:,44,k)=io(:,105,k)+io(:,110,k)+0.118*0.153*(io(:,108,k))+0.131*0.902*(io(:,109,k))+0.94*0.098*(io(:,109,k));
blstodj(:,45,k)=io(:,104,k)+0.165*(io(:,108,k));
blstodj(:,46,k)=io(:,102,k)+0.249*(io(:,108,k));
blstodj(:,47,k)=io(:,107,k)+io(:,141,k)+0.319*(io(:,108,k))+0.102*0.868*(io(:,109,k));
blstodj(:,48,k)=io(:,115,k)+0.017*(io(:,136,k));
blstodj(:,49,k)=io(:,114,k);
blstodj(:,50,k)=io(:,12,k)+io(:,178,k)+0.016*(io(:,136,k));
blstodj(:,51,k)=io(:,13,k);
blstodj(:,52,k)=io(:,14,k)+(io(:,145,k));
blstodj(:,53,k)=io(:,100,k)+(0.034+0.035)*(io(:,136,k));
blstodj(:,54,k)=io(:,101,k)-(0.122*0.94+0.061*0.68)*io(:,101,k)+(0.042+0.037+0.024+0.019+0.015)*(io(:,136,k));
blstodj(:,55,k)=0.013*io(:,110,k)+(0.122*0.94+0.061*0.68)*io(:,101,k)+(0.014)*(io(:,136,k));
blstodj(:,56,k)=io(:,164,k)+(0.034)*(io(:,136,k));
blstodj(:,57,k)=(0.009+(1-0.009)*0.719)*io(:,117,k)+0.028*blstodj(:,118,k)+(0.073)*(io(:,136,k));
blstodj(:,58,k)=io(:,126,k)+((1-0.009)*0.266)*io(:,117,k)+(0.915+0.04)*(io(:,118,k))+(0.459)*(io(:,121,k))+(0.444*0.338)*(io(:,174,k))+(0.01+0.063)*(io(:,136,k));
blstodj(:,59,k)=io(:,119,k)+io(:,120,k)+(0.015)*(io(:,118,k))+(0.541)*(io(:,121,k))+(0.033)*(io(:,136,k));
%blstodj(:,60,k)= real estate owner occupied imputed rents; not considered
blstodj(:,61,k)=0.977*io(:,122,k)+0.107*io(:,159,k)+(0.134*0.1)*(io(:,15,k))+0.05*sum(io(:,169:175,k)')'+(0.012)*(io(:,136,k));
%blstodj(:,62,k)= renting of machinery; no data in original jorgenson data
blstodj(:,63,k)=io(:,163,k);
blstodj(:,64,k)=0.825*sum(io(:,169:172,k)')'+0.027*sum(io(:,123:125,k)')'+0.03*sum(io(:,127:135,k)')';
blstodj(:,65,k)=0.286*sum(io(:,111:112,k)')'-io(:,112,k)+0.033*(io(:,113,k))+((1-0.674)*0.365)*(io(:,116,k))+0.382*sum(io(:,123:125,k)')'-0.088*(io(:,125,k))+0.284*sum(io(:,127:135,k)')'-0.962*(io(:,131,k))+0.837*sum(io(:,137:144,k)')'+0.02*(io(:,145,k))+0.012*sum(io(:,157:160,k)')'+0.023*sum(io(:,30:31,k)')'+0.058*sum(io(:,169:172,k)')'+0.027*io(:,42,k)+0.01*io(:,136,k);
blstodj(:,66,k)=io(:,112,k)+0.431*(io(:,83,k))+((0.674)*0.98)*(io(:,116,k))+0.032*0.04*io(:,101,k)+0.011*0.102*io(:,108,k)+0.088*(io(:,125,k))+0.962*(io(:,131,k))+0.03*sum(io(:,165:168,k)')'+0.015*io(:,136,k);
blstodj(:,67,k)=0.712*sum(io(:,165:168,k)')'+0.081*0.953*(io(:,108,k))+0.311*sum(io(:,123:125,k)')'+0.063*sum(io(:,169:172,k)')';
blstodj(:,68,k)=0.23*sum(io(:,165:168,k)')'+0.058*(io(:,145,k))+(0.011*0.036+0.026*0.038)*(io(:,101,k));
blstodj(:,69,k)=0.246*sum(io(:,123:125,k)')'+0.922*(io(:,113,k))+(0.018)*sum(io(:,157:160,k)')';
blstodj(:,70,k)=0.843*sum(io(:,157:160,k)')'+0.038*sum(io(:,146:148,k)')'+0.041*(io(:,161,k))+0.98*(io(:,162,k))+(0.376*0.047)*(io(:,108,k))+(0.022)*sum(io(:,123:125,k)')';
blstodj(:,71,k)=io(:,149,k)+(0.533*0.292)*(io(:,151,k))+(0.016)*(io(:,136,k));
blstodj(:,72,k)=0.582*io(:,153,k)+(0.171*0.15+0.186*0.443)*(io(:,153,k));
blstodj(:,73,k)=io(:,152,k)+(0.022)*(io(:,136,k));
blstodj(:,74,k)=io(:,150,k)+(0.15)*(io(:,98,k))+(0.226*0.989+0.241*0.426+0.708*0.533)*(io(:,151,k));;
blstodj(:,75,k)=0.944*io(:,127,k);
blstodj(:,76,k)=0.943*sum(io(:,146:148,k)')'+0.326*0.573*(io(:,116,k))+0.011*(io(:,136,k));
blstodj(:,77,k)=0.31*io(:,153,k)+sum(io(:,154:156,k)')'+(0.165+0.779)*sum(io(:,173:175,k)')'+0.959*(io(:,161,k))+(0.032+0.011)*(io(:,136,k));
blstodj(:,78,k)=0.084*io(:,129,k)+0.825*io(:,133,k)+(0.254)*(io(:,135,k));
blstodj(:,79,k)=0.445*sum(io(:,127:135,k)')'+(0.037)*sum(io(:,157:160,k)')'+0.057*0.326*io(:,116,k)+(0.048)*sum(io(:,137:144,k)')'+(0.01)*sum(io(:,146:148,k)')'+0.025*0.235*io(:,15,k)+0.014*0.134*io(:,15,k)-blstodj(:,78,k);
blstodj(:,80,k)=io(:,176,k);
blstodj(:,81,k)=sum(io(:,183:185,k)')';
blstodj(:,82,k)=io(:,177,k)+io(:,179,k);
blstodj(:,83,k)=io(:,188,k)+io(:,192,k);
%blstodj(:,84,k)= no match available
blstodj(:,85,k)=io(:,189,k)+io(:,193,k);
blstodj(:,86,k)=io(:,190,k)+sum(io(:,194:196,k)')';
blstodj(:,87,k)=io(:,191,k)+sum(io(:,186:187,k)')';
blstodj(:,88,k)=sum(io(:,180:182,k)')';

% Aggregating Rows of BLS IO data

blstodj(1,:,k)=sum(blstodj(1:2,:,k));
blstodj(2,:,k)=(blstodj(6,:,k));
blstodj(3,:,k)=(blstodj(3,:,k))+(blstodj(5,:,k));
blstodj(4,:,k)=0.16*sum(blstodj(8:10,:,k));
blstodj(5,:,k)=0.475*sum(blstodj(8:10,:,k));
blstodj(6,:,k)=0.336*sum(blstodj(8:10,:,k))+0.016*(blstodj(11,:,k));
blstodj(7,:,k)=(blstodj(7,:,k))+0.971*(blstodj(11,:,k));
%blstodj(8,:,k)=%blstodj(:,8,k)= uranium ore mining; no data in original
blstodj(9,:,k)=0.235*0.915*(blstodj(15,:,k))+0.109*(blstodj(145,:,k))+0.235*0.054*(blstodj(15,:,k))+0.134*0.83*(blstodj(15,:,k))+0.134*0.054*(blstodj(15,:,k))+0.632*0.984*(blstodj(15,:,k));
blstodj(10,:,k)=sum(blstodj(37:39,:,k))+(blstodj(4,:,k));
blstodj(11,:,k)=sum(blstodj(95:97,:,k));
blstodj(12,:,k)=sum(blstodj(53:56,:,k));
blstodj(13,:,k)=sum(blstodj(57:61,:,k));
blstodj(14,:,k)=sum(blstodj(62:70,:,k));
blstodj(15,:,k)=sum(blstodj(71:72,:,k))+sum(blstodj(74:77,:,k))+(0.81-0.159)*blstodj(73,:,k);
blstodj(16,:,k)=blstodj(78,:,k);
blstodj(17,:,k)=0.018*(blstodj(63,:,k))+0.975*(blstodj(84,:,k))+0.399*0.054*(blstodj(87,:,k))+0.023*(blstodj(90,:,k));
blstodj(18,:,k)=blstodj(80,:,k)+0.078*0.211*(blstodj(113,:,k));
blstodj(19,:,k)=sum(blstodj(85:87,:,k))+0.99*(blstodj(84,:,k));
blstodj(20,:,k)=(blstodj(79,:,k));
blstodj(21,:,k)=(blstodj(81,:,k))+(blstodj(83,:,k))+0.01*(blstodj(136,:,k));
blstodj(22,:,k)=sum(blstodj(88:90,:,k))+0.016*(blstodj(136,:,k));
blstodj(23,:,k)=(blstodj(91,:,k));
blstodj(24,:,k)=(blstodj(93,:,k));
blstodj(25,:,k)=(blstodj(92,:,k))+(blstodj(94,:,k));
blstodj(26,:,k)=blstodj(82,:,k)+0.19*(blstodj(73,:,k));
blstodj(27,:,k)=(blstodj(98,:,k));
blstodj(28,:,k)=0.314*(blstodj(50,:,k))+0.159*(blstodj(73,:,k))+0.01*(blstodj(84,:,k));
blstodj(29,:,k)=0.899*(blstodj(98,:,k));
blstodj(30,:,k)=sum(blstodj(16:25,:,k))+0.031*(blstodj(136,:,k));
blstodj(31,:,k)=(blstodj(26,:,k));
blstodj(32,:,k)=sum(blstodj(27:32,:,k));
blstodj(33,:,k)=0.982*(blstodj(33,:,k))+0.944*(blstodj(34,:,k));
blstodj(34,:,k)=0.056*(blstodj(34,:,k))+0.982*(blstodj(35,:,k))+0.89*(blstodj(36,:,k));
blstodj(35,:,k)=sum(blstodj(40:41,:,k))+0.013*(blstodj(136,:,k));
blstodj(36,:,k)=(blstodj(111,:,k))+0.078*0.126*(blstodj(113,:,k))+(1-0.699-0.07-0.057-0.048)*(blstodj(42,:,k));
blstodj(37,:,k)=(0.699+0.07+0.057+0.048)*(blstodj(42,:,k));
blstodj(38,:,k)=sum(blstodj(44:46,:,k))+sum(blstodj(48:49,:,k))+0.458*(blstodj(50,:,k))+0.014*(blstodj(136,:,k));
blstodj(39,:,k)=blstodj(47,:,k)+0.021*(blstodj(133,:,k));
blstodj(40,:,k)=(blstodj(43,:,k));
blstodj(41,:,k)=sum(blstodj(51:52,:,k))+0.11*(blstodj(36,:,k))+0.101*(blstodj(99,:,k));
blstodj(42,:,k)=blstodj(103,:,k)+0.102*0.11*(blstodj(108,:,k));
blstodj(43,:,k)=blstodj(106,:,k)+0.118*0.032*(blstodj(108,:,k));
blstodj(44,:,k)=blstodj(105,:,k)+blstodj(110,:,k)+0.118*0.153*(blstodj(108,:,k))+0.131*0.902*(blstodj(109,:,k))+0.94*0.098*(blstodj(109,:,k));
blstodj(45,:,k)=blstodj(104,:,k)+0.165*(blstodj(108,:,k));
blstodj(46,:,k)=blstodj(102,:,k)+0.249*(blstodj(108,:,k));
blstodj(47,:,k)=blstodj(107,:,k)+blstodj(141,:,k)+0.319*(blstodj(108,:,k))+0.102*0.868*(blstodj(109,:,k));
blstodj(48,:,k)=blstodj(115,:,k)+0.017*(blstodj(136,:,k));
blstodj(49,:,k)=blstodj(114,:,k);
blstodj(50,:,k)=blstodj(12,:,k)+blstodj(178,:,k)+0.016*(blstodj(136,:,k));
blstodj(51,:,k)=blstodj(13,:,k);
blstodj(52,:,k)=blstodj(14,:,k)+(blstodj(145,:,k));
blstodj(53,:,k)=blstodj(100,:,k)+(0.034+0.035)*(blstodj(136,:,k));
blstodj(54,:,k)=blstodj(101,:,k)-(0.122*0.94+0.061*0.68)*blstodj(101,:,k)+(0.042+0.037+0.024+0.019+0.015)*(blstodj(136,:,k));
blstodj(55,:,k)=0.013*blstodj(110,:,k)+(0.122*0.94+0.061*0.68)*blstodj(101,:,k)+(0.014)*(blstodj(136,:,k));
blstodj(56,:,k)=blstodj(164,:,k)+(0.034)*(blstodj(136,:,k));
blstodj(57,:,k)=(0.009+(1-0.009)*0.719)*blstodj(117,:,k)+0.028*blstodj(118,:,k)+(0.073)*(blstodj(136,:,k));
blstodj(58,:,k)=blstodj(126,:,k)+((1-0.009)*0.266)*blstodj(117,:,k)+(0.915+0.04)*(blstodj(118,:,k))+(0.459)*(blstodj(121,:,k))+(0.444*0.338)*(blstodj(174,:,k))+(0.01+0.063)*(blstodj(136,:,k));
blstodj(59,:,k)=blstodj(119,:,k)+blstodj(120,:,k)+(0.015)*(blstodj(118,:,k))+(0.541)*(blstodj(121,:,k))+(0.033)*(blstodj(136,:,k));
%blstodj(:,60,k)= real estate owner occupied; not considered...
blstodj(61,:,k)=0.977*blstodj(122,:,k)+0.107*blstodj(159,:,k)+(0.134*0.1)*(blstodj(15,:,k))+0.05*sum(blstodj(169:175,:,k))+(0.012)*(blstodj(136,:,k));
%blstodj(:,62,k)= renting of machinery; no data in original jorgenson data
blstodj(63,:,k)=blstodj(163,:,k);
blstodj(64,:,k)=0.825*sum(blstodj(169:172,:,k))+0.027*sum(blstodj(123:125,:,k))+0.03*sum(blstodj(127:135,:,k));
blstodj(65,:,k)=0.286*sum(blstodj(111:112,:,k))-blstodj(112,:,k)+0.033*(blstodj(113,:,k))+((1-0.674)*0.365)*(blstodj(116,:,k))+0.382*sum(blstodj(123:125,:,k))-0.088*(blstodj(125,:,k))+0.284*sum(blstodj(127:135,:,k))-0.962*(blstodj(131,:,k))+0.837*sum(blstodj(137:144,:,k))+0.02*(blstodj(145,:,k))+0.012*sum(blstodj(157:160,:,k))+0.023*sum(blstodj(30:31,:,k))+0.058*sum(blstodj(169:172,:,k))+0.027*blstodj(42,:,k)+0.01*blstodj(136,:,k);
blstodj(66,:,k)=blstodj(112,:,k)+(blstodj(83,:,k))+((0.674)*0.98)*(blstodj(116,:,k))+0.032*0.04*blstodj(101,:,k)+0.011*0.102*blstodj(108,:,k)+0.088*(blstodj(125,:,k))+0.962*(blstodj(131,:,k))+0.03*sum(blstodj(165:168,:,k))+0.015*blstodj(136,:,k);
blstodj(67,:,k)=0.712*sum(blstodj(165:168,:,k))+0.081*0.953*(blstodj(108,:,k))+0.311*sum(blstodj(123:125,:,k))+0.063*sum(blstodj(169:172,:,k));
blstodj(68,:,k)=0.23*sum(blstodj(165:168,:,k))+0.058*(blstodj(145,:,k))+(0.011*0.036+0.026*0.038)*(blstodj(101,:,k));
blstodj(69,:,k)=0.246*sum(blstodj(123:125,:,k))+0.922*(blstodj(113,:,k))+(0.018)*sum(blstodj(157:160,:,k));
blstodj(70,:,k)=0.843*sum(blstodj(157:160,:,k))+0.038*sum(blstodj(146:148,:,k))+0.041*(blstodj(161,:,k))+0.98*(blstodj(162,:,k))+(0.376*0.047)*(blstodj(108,:,k))+(0.022)*sum(blstodj(123:125,:,k));
blstodj(71,:,k)=blstodj(149,:,k)+(0.533*0.292)*(blstodj(151,:,k))+(0.016)*(blstodj(136,:,k));
blstodj(72,:,k)=0.582*blstodj(153,:,k)+(0.171*0.15+0.186*0.443)*(blstodj(153,:,k));
blstodj(73,:,k)=blstodj(152,:,k)+(0.022)*(blstodj(136,:,k));
blstodj(74,:,k)=blstodj(150,:,k)+(0.15)*(blstodj(98,:,k))+(0.226*0.989+0.241*0.426+0.708*0.533)*(blstodj(151,:,k));
blstodj(75,:,k)=0.944*blstodj(127,:,k);
blstodj(76,:,k)=0.943*sum(blstodj(146:148,:,k))+0.326*0.573*(blstodj(116,:,k))+0.011*(blstodj(136,:,k));
blstodj(77,:,k)=0.31*blstodj(153,:,k)+sum(blstodj(154:156,:,k))+(0.165+0.779)*sum(blstodj(173:175,:,k))+0.959*(blstodj(161,:,k))+(0.032+0.011)*(blstodj(136,:,k));
blstodj(78,:,k)=0.084*blstodj(129,:,k)+0.825*blstodj(133,:,k)+(0.254)*(blstodj(135,:,k));
blstodj(79,:,k)=0.445*sum(blstodj(127:135,:,k))+(0.037)*sum(blstodj(157:160,:,k))+0.057*0.326*blstodj(116,:,k)+(0.048)*sum(blstodj(137:144,:,k))+(0.01)*sum(blstodj(146:148,:,k))+0.025*0.235*blstodj(15,:,k)+0.014*0.134*blstodj(15,:,k)-blstodj(78,:,k);
blstodj(80,:,k)=blstodj(176,:,k);
blstodj(81,:,k)=sum(blstodj(183:185,:,k));
blstodj(82,:,k)=blstodj(177,:,k)+blstodj(179,:,k);
blstodj(83,:,k)=blstodj(188,:,k)+blstodj(192,:,k);
%blstodj(:,84,k)= no match available
blstodj(85,:,k)=blstodj(189,:,k)+blstodj(193,:,k);
blstodj(86,:,k)=blstodj(190,:,k)+sum(blstodj(194:196,:,k));
blstodj(87,:,k)=blstodj(191,:,k)+sum(blstodj(186:187,:,k));
blstodj(88,:,k)=sum(blstodj(180:182,:,k));
end

%88 Sector x Time Panel: Aggregated from BLS data

ioblstodj=[blstodj(1:88,1:88,:);blstodj(203,1:88,:)];

for i=1:size(ioblstodj,3)
    for j=1:size(ioblstodj,2)
        vaddblsdja(j,i)=(ioblstodj(end,j,i));
    end
end

for i=1:size(ioblstodj,3)
    for j=1:size(ioblstodj,2)
        grossyblsdja(j,i)=sum(ioblstodj(:,j,i));
    end
end

%REMOVE GOVERNMENT SECTORS & RENTS IMPUTED FROM OWNER-OCCUPIED HOUSING
temp=81:1:88;
temp=[60;temp'];
grossyblsdja(temp,:)=[];
vaddblsdja(temp,:)=[];
clear temp

%REMOVE SECTORS FOR WHICH THERE IS NO GROSS SALES (8, URANIUM ORES AND 62,
%RENTING OF MACHINERY)

[temp]=find(sum(grossyblsdja')==0);
grossyblsdja(temp,:)=[];
vaddblsdja(temp,:)=[];

gdpblstodja=sum(vaddblsdja);

%% 3) SPLICE BLS DATA WITH JORGENSON SECTORS USING GROWTH RATES IN BLS DATA

for i=1:size(grossyblsdja,1)
    for j=2:size(grossyblsdja,2)
        growthgobls(i,j)=(grossyblsdja(i,j)-grossyblsdja(i,j-1))/grossyblsdja(i,j-1);
        if isnan(growthgobls(i,j))==1
            growthgobls(i,j)=0;
        end
    end
end

for i=1:size(vaddblsdja,1)
    for j=2:size(vaddblsdja,2)
        growthvaddbls(i,j)=(vaddblsdja(i,j)-vaddblsdja(i,j-1))/vaddblsdja(i,j-1);
        if isnan(growthvaddbls(i,j))==1
            growthvaddbls(i,j)=0
        end
    end
end


% SPLICE USING 1999 as BASE YEAR
% 1999=40 in old dja; 7 in new bls
grossalt=[grossy(:,1:40) zeros(size(grossy,1),(2008-1999))];
for i=1:size(grossalt,1)
    for j=1:(2008-1999)
        grossalt(i,40+j)=(grossalt(i,40+j-1))*(1+growthgobls(i,7+j));
    end
end

vaddalt=[vadd(:,1:40) zeros(size(vadd,1),(2008-1999))];
for i=1:size(vaddalt,1)
    for j=1:(2008-1999)
        vaddalt(i,40+j)=(vaddalt(i,40+j-1))*(1+growthvaddbls(i,7+j));
    end
end

 %Compute new Domar Weights
 for i=1:size(vaddalt,1)
     for j=1:size(vaddalt,2)
         dweightalt(i,j)=grossalt(i,j)/(sum(vaddalt(:,j)));
     end
 end

hdalt=sum(dweightalt.^2)


