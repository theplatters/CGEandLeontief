%% US DATA SERIES FOR THE GREAT DIVERSIFICATION AND ITS UNDOING
% PRIVATE SECTOR DATA 1960 - 2008
%1) Takes Original 88 Industries by Jorgenson et al. and Removes Government Sectors + Imputed
%rents to owner occupied households + sectors for which there is no data in
%original (Uranium Ore and Machine Rental)
%2) Computes TFP accounting for all sectors and Aggregates
%3) Splices Jorgenson et al data (1960-2005) with BLS data till 2008
%% 

%% 1) READ JORGENSON ET AL DATA AND PERFORM BASIC MANIPULATIONS
%data is quantity data matrix
%price is price data matrix
clear
load us80dbasedata.mat

%AUXILIARY COUNTERS TO TURN ORIGINAL DATA INTO MATRIX(SECTORxYEAR)
startcount=1:46:4048;
for k=1:size(startcount,2)
    endcount(k)=startcount(k)-1;
end
endcount(89)=size(data,1);

%NOMINAL GROSS OUTPUT IS THE SUM OF NOMINAL CAPITAL, LABOR AND ALL
%INTERMEDIATE INPUTS (INCLUDING NONCOMPETING IMPORTS)
for i=1:size(startcount,2)
    grossy(i,:)=data(startcount(i):endcount(i+1),3)';
end

%NOMINAL VALUE ADDED IS THE SUM OF NOMINAL CAPITAL AND LABOR
for i=1:size(startcount,2)
    vadd(i,:)=(data(startcount(i):endcount(i+1),4)')+(data(startcount(i):endcount(i+1),5)');
end

%NOMINAL CAPITAL
for i=1:size(startcount,2)
    capital(i,:)=(data(startcount(i):endcount(i+1),4)');
end

%NOMINAL LABOR
for i=1:size(startcount,2)
    labor(i,:)=(data(startcount(i):endcount(i+1),5)');
end

%REMOVE GOVERNMENT SECTORS & RENTS IMPUTED FROM OWNER-OCCUPIED HOUSING
temp=81:1:88;
temp=[60;temp'];
grossy(temp,:)=[];
vadd(temp,:)=[];
capital(temp,:)=[];
labor(temp,:)=[];
clear temp


%REMOVE SECTORS FOR WHICH THERE IS NO GROSS SALES (8, URANIUM ORES AND 62,
%RENTING OF MACHINERY)

[temp]=find(sum(grossy')==0);
grossy(temp,:)=[];
vadd(temp,:)=[];
capital(temp,:)=[];
labor(temp,:)=[];

%VALUE ADDED WEIGHTS IS VA(I,J)/SUM(VA(:,J))
for i=1:size(vadd,1)
   for j=1:size(vadd,2)
       vaddshare(i,j)=vadd(i,j)/(sum(vadd(:,j)));
   end
end 

%DOMAR WEIGHTS IS GO(I,J)/SUM(VA(:,J))
for i=1:size(vadd,1)
    for j=1:size(vadd,2)
        dweight(i,j)=grossy(i,j)/(sum(vadd(:,j)));
    end
end

%CAPITAL WEIGHT IN GROSS OUTPUT
for i=1:size(vadd,1)
    for j=1:size(vadd,2)
        kweight(i,j)=capital(i,j)/((grossy(i,j)));
    end
end

%LABOR WEIGHT IN GROSS OUTPUT
for i=1:size(vadd,1)
    for j=1:size(vadd,2)
        lweight(i,j)=labor(i,j)/((grossy(i,j)));
    end
end


%% 2) TFP ACCOUNTING (EQT 7.08 P296 IN JORGENSON ET AL 2005)
%VOLUME INDICES first
%ALL PRICES ARE 1996=1

vol=data./price;
for i=2:size(vol,1)
    dvol(i,:)=log(vol(i,:))-log(vol(i-1,:));
end

%REAL GROSS OUTPUT GROWTH
for i=1:size(startcount,2)
    dvoly(i,:)=dvol(startcount(i):endcount(i+1),3)';
end

%REAL CAPITAL STOCK GROWTH
for i=1:size(startcount,2)
    dvolk(i,:)=dvol(startcount(i):endcount(i+1),4)';
end

%REAL LABOUR GROWTH
for i=1:size(startcount,2)
    dvoll(i,:)=dvol(startcount(i):endcount(i+1),5)';
end

%NOTE:zero out first columns for 1960
dvoly(:,1)=0;
dvolk(:,1)=0;
dvoll(:,1)=0;

%REAL INTERMEDIATE INPUT GROWTH
%sequence of matrices (input by year by sector)
%number of inputs, including noncompetitive imports is size(dvol,2)
for i=1:size(startcount,2)
    dvolii(:,:,i)=dvol(startcount(i):endcount(i+1),6:size(dvol,2))';
end

%NOTE:zero out first columns for 1960
dvolii(:,1,:)=0;

%REMOVE GOVERNMENT SECTORS
clear temp
temp=81:1:88;
temp=[60;temp'];
dvoly(temp,:)=[];
dvolk(temp,:)=[];
dvoll(temp,:)=[];
dvolii(:,:,temp)=[];
clear temp


%REMOVE SECTORS FOR WHICH THERE IS NO GROSS SALES (8, URANIUM ORES AND 62,
%RENTING OF MACHINERY)

[temp]=[8;61];
dvoly(temp,:)=[];
dvolk(temp,:)=[];
dvoll(temp,:)=[];
dvolii(:,:,temp)=[];

%TWO PERIOD AVERAGE VALUE SHARES IN NOMINAL WEIGHTS
%TWO PERIOD AVERAGE CAPITAL and LABOR SHARES IN GROSS OUTPUT
for i=1:size(vadd,1)
   for j=2:size(vadd,2)
       ksharebar(i,j)=0.5*(kweight(i,j)+kweight(i,j-1));
       lsharebar(i,j)=0.5*(lweight(i,j)+lweight(i,j-1));
   end
end 

%TWO PERIOD AVERAGE INTERMEDIATE INPUT SHARES IN GROSS OUTPUT
%FIRST FIND WEIGHTS FOR EACH INTERMEDIATE IN TOTAL INTERMEDIATE INPUT BILL
%OF EACH SECTOR

%TOTAL INTERMEDIATE INPUT BILL
totii=grossy-vadd;
%very small negative values due to discrepancies (set them to zero)
for i=1:size(totii,1)
    for j=1:size(totii,2)
        if totii(i,j)<1
            totii(i,j)=0;
        end
    end
end

%INTERMEDIATE INPUT MATRICES BY YEAR BY SECTOR
%sequence of matrices (input by year by sector)
%number of inputs, including noncompetitive imports is size(dvol,2)
for i=1:size(startcount,2)
    ii(:,:,i)=data(startcount(i):endcount(i+1),6:size(data,2))';
end


clear temp
temp=81:1:88;
temp=[60;temp'];
ii(:,:,temp)=[];
clear temp

[temp]=[8;61];
ii(:,:,temp)=[];

%TWO PERIOD AVERAGE WEIGHTS

for i=1:size(ii,3)
    for j=1:size(ii,2)
        for k=1:size(ii,1)
            iiweight(k,j,i)=ii(k,j,i)/totii(i,j);
        end
    end
end

for i=1:size(ii,3)
    for j=2:size(ii,2)
        for k=1:size(ii,1)
            iisharebar(k,j,i)=0.5*(iiweight(k,j,i)+iiweight(k,j-1,i));
        end
    end
end

for i=1:size(iisharebar,1)
    for j=1:size(iisharebar,2)
        for k=1:size(iisharebar,3)
            if isnan(iisharebar(i,j,k))==1
               iisharebar(i,j,k)=0;
            end
        end
    end
end

for i=1:size(dvolii,1)
    for j=1:size(dvolii,2)
        for k=1:size(dvolii,3)
            if isnan(dvolii(i,j,k))==1
                dvolii(i,j,k)=0;
            end
        end
    end
end

for i=1:size(dvolii,1)
    for j=1:size(dvolii,2)
        for k=1:size(dvolii,3)
            if isinf(dvolii(i,j,k))==1
                dvolii(i,j,k)=0;
            end
        end
    end
end

%TORNQVIST INDEX OF INTERMEDIATE INPUT WEIGHTS
%See Jorgenson et al (2005) pg 97, eqn 4.04

for i=1:size(ii,3)
    for j=2:size(ii,2)
        dvoltotii(i,j)=(iisharebar(:,j,i)')*(dvolii(:,j,i));
    end
end

%1996 is base year (column 37 of dvoltotii)
%Backing out price and real quantity of intermediate input as in Jorgenson
%et al (2005), eqn 4.07

xi(:,37)=totii(:,37);

for i=1:size(dvoltotii,1)
    for j=38:size(dvoltotii,2)
        xi(i,j)=xi(i,j-1)*(exp(dvoltotii(i,j)));
    end
end

for i=1:size(dvoltotii,1)
    for j=1:36
        xi(i,37-j)=xi(i,37-j+1)/(exp(dvoltotii(i,37-j+1)));
    end
end

for i=1:size(dvoltotii,1)
    for j=1:size(dvoltotii,2)
        pii(i,j)=totii(i,j)/xi(i,j);
    end
end

%VALUE SHARE OF INTERMEDIATE INPUTS

for i=1:size(dvoltotii,1)
    for j=1:size(dvoltotii,2)
        totiishare(i,j)=(pii(i,j)*xi(i,j))/grossy(i,j);
    end
end

for i=1:size(dvoltotii,1)
    for j=2:size(dvoltotii,2)
        totiisharebar(i,j)=0.5*(totiishare(i,j)+totiishare(i,j-1));
    end
end

%SECTORAL TFP GROWTH (Gross Output Perspective)
%FORMULA 7.08 in Jorgenson et al (2005) page 296 
for i=1:size(dvoltotii,1)
    for j=2:size(dvoltotii,2)
        stfp(i,j)=dvoly(i,j)-ksharebar(i,j)*dvolk(i,j)-lsharebar(i,j)*dvoll(i,j)-totiisharebar(i,j)*dvoltotii(i,j);
    end
end

for i=1:size(stfp,1)
    for j=1:size(stfp,2)
        if isnan(stfp(i,j))==1
            stfp(i,j)=0;
        end
    end
end

%SECTORAL TFP GROWTH (Value Added Perspective)
%FORMULA 7.16 in Jorgenson et al (2005) page 296 

%First need two-period average of value added weight in gross output each sector

for i=1:size(vadd,1)
   for j=1:size(vadd,2)
       vaddweight(i,j)=vadd(i,j)/(grossy(i,j));
   end
end 

for i=1:size(vaddweight,1)
    for j=2:size(vaddweight,2)
        vaddweightbar(i,j)=0.5*(vaddweight(i,j)+vaddweight(i,j-1));
    end
end

%Now FORMULA 7.16 in Jorgenson et al (2005) page 296 
for i=1:size(stfp,1)
    for j=2:size(stfp,2)
        stfpva(i,j)=stfp(i,j)/vaddweightbar(i,j);
        stfpva2(i,j)=stfp(i,j)/vaddweight(i,j);
    end
end



%AGGREGATE TFP GROWTH
%WITH 2 PERIOD AVERAGE DOMAR WEIGHTS

for i=1:size(dweight,1)
    for j=2:size(dweight,2)
        dweightbar(i,j)=0.5*(dweight(i,j)+dweight(i,j-1));
    end
end

for j=2:size(dweight,2)
        aggtfp(j)=(stfp(:,j)'*dweightbar(:,j));
end


%AGGREGATE TFP STATISTICS
meanaggtfp=mean(aggtfp(2:end).*100);
varaggtfp=var(aggtfp(2:end).*100);
prevaraggtfp=var(aggtfp(2:23).*100);
postvaraggtfp=var(aggtfp(24:end).*100);
postpreratio=postvaraggtfp/prevaraggtfp;
stdaggtfp=std(aggtfp(2:end).*100);
prestdaggtfp=std(aggtfp(2:23).*100);
poststdaggtfp=std(aggtfp(24:end).*100);
postpreratiostd=poststdaggtfp/prestdaggtfp;


%AGGREGATE OUTPUT GROWTH
%Value added shares, 2 period averages
for i=1:size(vaddshare,1)
    for j=2:size(vaddshare,2)
        vaddsharebar(i,j)=0.5*(vaddshare(i,j)+vaddshare(i,j-1));
    end
end

%Share of industry value added in industry gross-output
for i=1:size(vadd,1)
    for j=1:size(vadd,2)
        vweight(i,j)=vadd(i,j)/((grossy(i,j)));
    end
end
%2 period average of Share of industry value added in industry gross-output
for i=1:size(vweight,1)
    for j=2:size(vweight,2)
        dvweightbar(i,j)=0.5*(vweight(i,j)+vweight(i,j-1));
    end
end

for i=1:size(dvolk,1)
    for j=1:size(dvolk,2)
            if isnan(dvolk(i,j))==1
                dvolk(i,j)=0;
            end
    end
end

%Aggregate GDP GROWTH (formula 8.33 in Jorgenson et al (2005))
for j=2:size(dweight,2)
    agggdp(j)=0;
        for i=1:size(dweight,1)
            agggdp(j)=vaddsharebar(i,j)*ksharebar(i,j)*(1/dvweightbar(i,j))*dvolk(i,j)+vaddsharebar(i,j)*lsharebar(i,j)*(1/dvweightbar(i,j))*dvoll(i,j)+vaddsharebar(i,j)*(1/dvweightbar(i,j))*stfp(i,j)+agggdp(j);
        end
end

%AGGREGATE GDP STATISTICS
meanagggdp=mean(agggdp(2:end).*100);
varagggdp=var(agggdp(2:end).*100);
prevaragggdp=var(agggdp(2:23).*100);
postvaragggdp=var(agggdp(24:end).*100);
postpreratio=postvaragggdp/prevaragggdp;
stdagggdp=std(agggdp(2:end).*100);
prestdagggdp=std(agggdp(2:25).*100);
poststdagggdp=std(agggdp(26:end).*100);
postpreratiostd=poststdagggdp/prestdagggdp;

%% SPLICE WITH BLS DATA TO OBTAIN 2006-2008 DATA

[grossalt, vaddalt, dweightalt]=BLS_to_Jorgenson(grossy, vadd);

%save 'allUSdata.mat;