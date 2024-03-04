%% Two dimensional
% clear;
clc;
% Az in be baad be jaye inke x, y, z tarif koni 
% mituni az "imread" estefade koni
% x=xlsread('blah','a:a');
% y=xlsread('blah','b:b');
% z=xlsread('blah','c:c');
% [X,Y,Z] = xyz2grid(x,y,z);
Z=double(imread('KDGX_N0R_20130130_103700.tif'));
% img=imagesc(Z);
%% http://wradlib.org/wradlib-docs/0.6.0/tutorial_get_rainfall.html
% now we have our Reflectivity(Z) in logarithmic sApplevelle of 10 so its time to
% change this into rainfall intensity via various methods
% An Integrated Approach to Error Correction for Real-Time Radar-Rainfall
% Estimation for 53 and 15
% Or "Effects of underrepresented hydrometeor variability and partial beam
% filling on microwave brightness temperatures for rainfall retrieval" for
% 53 and 20 tu male Harris ham has
% Z(Z>53)=53;
% Z(Z<20)=0;
% Z(isnan(Z))=0;
% Z=Z(1:480,1:784);
%% in ghesmat baraye dorost kardan size anjam shode
x=ceil(size(Z,1)/128)*128-size(Z,1);
y=ceil(size(Z,2)/128)*128-size(Z,2);
left=ones(size(Z,1),floor(y/2))*min(min(Z));
right=ones(size(Z,1),ceil(y/2))*min(min(Z));
up=ones(floor(x/2),ceil(size(Z,2)/128)*128)*min(min(Z));
down=ones(ceil(x/2),ceil(size(Z,2)/128)*128)*min(min(Z));

Z=[left Z right];
Z=[up;Z;down];

a=200;
b=1.6;
R=nthroot((10.^(Z/10))/a,b);
R(R<0.5)=0;

% [a1,L1]=wavedec2(R,1,'db1');
% [a2,L2]=wavedec2(R,2,'db1');
% [a3,L3]=wavedec2(R,3,'db1');
% [a4,L4]=wavedec2(R,4,'db1');

[cA1,cH1,cV1,cD1]=dwt2(R,'haar');
[cA2,cH2,cV2,cD2]=dwt2(cA1,'haar');
[cA3,cH3,cV3,cD3]=dwt2(cA2,'haar');
[cA4,cH4,cV4,cD4]=dwt2(cA3,'haar');
[cA5,cH5,cV5,cD5]=dwt2(cA4,'haar');
[cA6,cH6,cV6,cD6]=dwt2(cA5,'haar');
[cA7,cH7,cV7,cD7]=dwt2(cA6,'haar');
%% in order to exclude non-related zeros 
for i=1:size(cA1,1)
    for j=1:size(cA1,2)
        if [cA1(i,j),cH1(i,j),cV1(i,j),cD1(i,j)]==[0,0,0,0]
            cA1(i,j)=NaN;
            cD1(i,j)=NaN;
            cH1(i,j)=NaN;
            cV1(i,j)=NaN;
        end
    end
end

for i=1:size(cA2,1)
    for j=1:size(cA2,2)
        if [cA2(i,j),cH2(i,j),cV2(i,j),cD2(i,j)]==[0,0,0,0]
            cA2(i,j)=NaN;
            cD2(i,j)=NaN;
            cH2(i,j)=NaN;
            cV2(i,j)=NaN;
        end
    end
end

for i=1:size(cA3,1)
    for j=1:size(cA3,2)
        if [cA3(i,j),cH3(i,j),cV3(i,j),cD3(i,j)]==[0,0,0,0]
            cA3(i,j)=NaN;
            cD3(i,j)=NaN;
            cH3(i,j)=NaN;
            cV3(i,j)=NaN;
        end
    end
end

for i=1:size(cA4,1)
    for j=1:size(cA4,2)
        if [cA4(i,j),cH4(i,j),cV4(i,j),cD4(i,j)]==[0,0,0,0]
            cA4(i,j)=NaN;
            cD4(i,j)=NaN;
            cH4(i,j)=NaN;
            cV4(i,j)=NaN;
        end
    end
end
for i=1:size(cA5,1)
    for j=1:size(cA5,2)
        if [cA5(i,j),cH5(i,j),cV5(i,j),cD5(i,j)]==[0,0,0,0]
            cA5(i,j)=NaN;
            cD5(i,j)=NaN;
            cH5(i,j)=NaN;
            cV5(i,j)=NaN;
        end
    end
end
for i=1:size(cA6,1)
    for j=1:size(cA6,2)
        if [cA6(i,j),cH6(i,j),cV6(i,j),cD6(i,j)]==[0,0,0,0]
            cA6(i,j)=NaN;
            cD6(i,j)=NaN;
            cH6(i,j)=NaN;
            cV6(i,j)=NaN;
        end
    end
end
for i=1:size(cA7,1)
    for j=1:size(cA7,2)
        if [cA7(i,j),cH7(i,j),cV7(i,j),cD7(i,j)]==[0,0,0,0]
            cA7(i,j)=NaN;
            cD7(i,j)=NaN;
            cH7(i,j)=NaN;
            cV7(i,j)=NaN;
        end
    end
end
epsilon5diag=Divide(cA5,cD5);
epsilon6diag=Divide(cA6,cD6);
epsilon7diag=Divide(cA7,cD7);
epsilon4diag=Divide(cA4,cD4);
epsilon3diag=Divide(cA3,cD3);
epsilon2diag=Divide(cA2,cD2);
epsilon1diag=Divide(cA1,cD1);

epsilon5hor=Divide(cA5,cH5);
epsilon6hor=Divide(cA6,cH6);
epsilon7hor=Divide(cA7,cH7);
epsilon4hor=Divide(cA4,cH4);
epsilon3hor=Divide(cA3,cH3);
epsilon2hor=Divide(cA2,cH2);
epsilon1hor=Divide(cA1,cH1);

epsilon5vert=Divide(cA5,cV5);
epsilon6vert=Divide(cA6,cV6);
epsilon7vert=Divide(cA7,cV7);
epsilon4vert=Divide(cA4,cV4);
epsilon3vert=Divide(cA3,cV3);
epsilon2vert=Divide(cA2,cV2);
epsilon1vert=Divide(cA1,cV1);

% stable1diag=fitdist(epsilon1diag(:),'stable');
% stable2diag=fitdist(epsilon2diag(:),'stable');
% stable3diag=fitdist(epsilon3diag(:),'stable');
% stable4diag=fitdist(epsilon4diag(:),'stable');

% Alpha_diag=[stable1diag.alpha,stable2diag.alpha,stable3diag.alpha,stable4diag.alpha];
% Beta_diag=[stable1diag.beta,stable2diag.beta,stable3diag.beta,stable4diag.beta];
% Gam_diag=[stable1diag.gam,stable2diag.gam,stable3diag.gam,stable4diag.gam];

% stable1hor=fitdist(epsilon1hor(:),'stable');
% stable2hor=fitdist(epsilon2hor(:),'stable');
% stable3hor=fitdist(epsilon3hor(:),'stable');
% stable4hor=fitdist(epsilon4hor(:),'stable');
% 
% Alpha_hor=[stable1hor.alpha,stable2hor.alpha,stable3hor.alpha,stable4hor.alpha];
% Beta_hor=[stable1hor.beta,stable2hor.beta,stable3hor.beta,stable4hor.beta];
% Gam_hor=[stable1hor.gam,stable2hor.gam,stable3hor.gam,stable4hor.gam];
% 
% stable1vert=fitdist(epsilon1vert(:),'stable');
% stable2vert=fitdist(epsilon3vert(:),'stable');
% stable3vert=fitdist(epsilon4vert(:),'stable');
% stable4vert=fitdist(epsilon4vert(:),'stable');
% 
% Alpha_vert=[stable1vert.alpha,stable2vert.alpha,stable3vert.alpha,stable4vert.alpha];
% Beta_vert=[stable1vert.beta,stable2vert.beta,stable3vert.beta,stable4vert.beta];
% Gam_vert=[stable1vert.gam,stable2vert.gam,stable3vert.gam,stable4vert.gam];
% tuye 1dimension az commande "fitdist" estefade mikonim ta inke enheraf meyar ha ro peyda
% konim vali inja chon nemishe fit kard be ellate matrix budanesh pas ba
% estefade az (:) tabdil be vector mikonim
% Semilogy behtarin halat va3 fit kardan hastesh
% 1-sum(sum(R))/sum(sum(r))
Fit5diag=fitdist(epsilon5diag(:),'normal');
% Fit6diag=fitdist(epsilon6diag(:),'normal');
% Fit7diag=fitdist(epsilon7diag(:),'normal');
Fit4diag=fitdist(epsilon4diag(:),'normal');
Fit3diag=fitdist(epsilon3diag(:),'normal');
Fit2diag=fitdist(epsilon2diag(:),'normal');
Fit1diag=fitdist(epsilon1diag(:),'normal');

Fit5vert=fitdist(epsilon5vert(:),'normal');
% Fit6vert=fitdist(epsilon6vert(:),'normal');
% Fit7vert=fitdist(epsilon7vert(:),'normal');
Fit4vert=fitdist(epsilon4vert(:),'normal');
Fit3vert=fitdist(epsilon3vert(:),'normal');
Fit2vert=fitdist(epsilon2vert(:),'normal');
Fit1vert=fitdist(epsilon1vert(:),'normal');

Fit5hor=fitdist(epsilon5hor(:),'normal');
% Fit6hor=fitdist(epsilon6hor(:),'normal');
% Fit7hor=fitdist(epsilon7hor(:),'normal');
Fit4hor=fitdist(epsilon4hor(:),'normal');
Fit3hor=fitdist(epsilon3hor(:),'normal');
Fit2hor=fitdist(epsilon2hor(:),'normal');
Fit1hor=fitdist(epsilon1hor(:),'normal');

Diag=[Fit1diag.sigma Fit2diag.sigma Fit3diag.sigma Fit4diag.sigma ... 
    Fit5diag.sigma];
Hor=[Fit1hor.sigma Fit2hor.sigma Fit3hor.sigma Fit4hor.sigma ... 
    Fit5hor.sigma];
Vert=[Fit1vert.sigma Fit2vert.sigma Fit3vert.sigma Fit4vert.sigma ...
    Fit5vert.sigma];

%% FIX THE FUCKING CODE instead of a for create fit use lowest resolution STD
% 3ta create fit function besaz "a" haro taghir bede be Fit1...sigma ha
% Sigma1Diag=Fit1diag.sigma;
% Sigma1Hor=Fit1hor.sigma;
% Sigma1Vert=Fit1vert.sigma;
m=[1,2,3,4,5];
[fitresultdiag, gof1] = CreateFitDiag(m, Diag);
[fitresulthor, gof2] = CreateFitHor(m, Hor);
[fitresultvert, gof3] = CreateFitVert(m, Vert);
%% Generating FluAppleveltuations
% STD4diag=fitresultdiag.a*2^(4*fitresultdiag.H);
% STD3diag=fitresultdiag.a*2^(3*fitresultdiag.H);
% STD2diag=fitresultdiag.a*2^(2*fitresultdiag.H);
% STD1diag=fitresultdiag.a*2^(1*fitresultdiag.H);
% 
% STD4hor=fitresulthor.a*2^(4*fitresulthor.H);
% STD3hor=fitresulthor.a*2^(3*fitresulthor.H);
% STD2hor=fitresulthor.a*2^(2*fitresulthor.H);
% STD1hor=fitresulthor.a*2^(1*fitresulthor.H);
% 
% STD4vert=fitresultvert.a*2^(4*fitresultvert.H);
% STD3vert=fitresultvert.a*2^(3*fitresultvert.H);
% STD2vert=fitresultvert.a*2^(2*fitresultvert.H);
% STD1vert=fitresultvert.a*2^(1*fitresultvert.H);
STD5diag=Fit5diag.sigma;
STD5hor=Fit5hor.sigma;
STD5vert=Fit5vert.sigma;

STD6diag=STD7diag*2^(-fitresultdiag.H);
STD6hor=STD7diag*2^(-fitresulthor.H);
STD6vert=STD7diag*2^(-fitresultvert.H);

STD5diag=STD7diag*2^(-fitresultdiag.H*2);
STD5hor=STD7diag*2^(-fitresulthor.H*2);
STD5vert=STD7diag*2^(-fitresultvert.H*2);

STD4diag=STD7diag*2^(-fitresultdiag.H*3);
STD4hor=STD7diag*2^(-fitresulthor.H*3);
STD4vert=STD7diag*2^(-fitresultvert.H*3);

STD3diag=STD7diag*2^(-fitresultdiag.H*4);
STD3hor=STD7diag*2^(-fitresulthor.H*4);
STD3vert=STD7diag*2^(-fitresultvert.H*4);

STD2diag=STD7diag*2^(-fitresultdiag.H*5);
STD2hor=STD7diag*2^(-fitresulthor.H*5);
STD2vert=STD7diag*2^(-fitresultvert.H*5);

STD1diag=STD7diag*2^(-fitresultdiag.H*6);
STD1hor=STD7diag*2^(-fitresulthor.H*6);
STD1vert=STD7diag*2^(-fitresultvert.H*6);

genstan1diag=normrnd(0,STD1diag,size(cA1));
genstan1diag(genstan1diag<-1)=-1;
genstan1diag(genstan1diag>1)=1;

genstan1hor=normrnd(0,STD2hor,size(cA1));
genstan1hor(genstan1hor<-1)=-1;
genstan1hor(genstan1hor>1)=1;

genstan1vert=normrnd(0,STD2vert,size(cA1));
genstan1vert(genstan1vert<-1)=-1;
genstan1vert(genstan1vert>1)=1;

genstan2diag=normrnd(0,STD2diag,size(cA2));
genstan2diag(genstan2diag<-1)=-1;
genstan2diag(genstan2diag>1)=1;

genstan2hor=normrnd(0,STD2hor,size(cA2));
genstan2hor(genstan2hor<-1)=-1;
genstan2hor(genstan2hor>1)=1;

genstan2vert=normrnd(0,STD2vert,size(cA2));
genstan2vert(genstan2vert<-1)=-1;
genstan2vert(genstan2vert>1)=1;

genstan3diag=normrnd(0,STD3diag,size(cA3));
genstan3diag(genstan3diag<-1)=-1;
genstan3diag(genstan3diag>1)=1;

genstan3hor=normrnd(0,STD3hor,size(cA3));
genstan3hor(genstan3hor<-1)=-1;
genstan3hor(genstan3hor>1)=1;

genstan3vert=normrnd(0,STD3vert,size(cA3));
genstan3vert(genstan3vert<-1)=-1;
genstan3vert(genstan3vert>1)=1;

genstan4diag=normrnd(0,STD4diag,size(cA4));
genstan4diag(genstan4diag<-1)=-1;
genstan4diag(genstan4diag>1)=1;

genstan4hor=normrnd(0,STD4hor,size(cA4));
genstan4hor(genstan4hor<-1)=-1;
genstan4hor(genstan4hor>1)=1;

genstan4vert=normrnd(0,STD4vert,size(cA4));
genstan4vert(genstan4vert<-1)=-1;
genstan4vert(genstan4vert>1)=1;

genstan5diag=normrnd(0,STD5diag,size(cA5));
genstan5diag(genstan5diag<-1)=-1;
genstan5diag(genstan5diag>1)=1;

genstan5hor=normrnd(0,STD5hor,size(cA5));
genstan5hor(genstan5hor<-1)=-1;
genstan5hor(genstan5hor>1)=1;

genstan5vert=normrnd(0,STD5vert,size(cA5));
genstan5vert(genstan5vert<-1)=-1;
genstan5vert(genstan5vert>1)=1;

genstan6diag=normrnd(0,STD6diag,size(cA6));
genstan6diag(genstan6diag<-1)=-1;
genstan6diag(genstan6diag>1)=1;

genstan6hor=normrnd(0,STD6hor,size(cA6));
genstan6hor(genstan6hor<-1)=-1;
genstan6hor(genstan6hor>1)=1;

genstan6vert=normrnd(0,STD6vert,size(cA6));
genstan6vert(genstan6vert<-1)=-1;
genstan6vert(genstan6vert>1)=1;

genstan7diag=normrnd(0,STD7diag,size(cA7));
genstan7diag(genstan7diag<-1)=-1;
genstan7diag(genstan7diag>1)=1;

genstan7hor=normrnd(0,STD7hor,size(cA7));
genstan7hor(genstan7hor<-1)=-1;
genstan7hor(genstan7hor>1)=1;

genstan7vert=normrnd(0,STD7vert,size(cA7));
genstan7vert(genstan7vert<-1)=-1;
genstan7vert(genstan7vert>1)=1;
%% now we start the downsApplevelling process
gendet7vert=genstan7vert.*cA7;
gendet7hor=genstan7hor.*cA7;
gendet7diag=genstan7diag.*cA7;
Appgen6=idwt2(cA7,gendet7hor,gendet7vert,gendet7diag,'haar');
% Appgen3(isnan(Appgen3))=0;
Appgen6=RemoveNeg(Appgen6);
% Appgen3=ArrangeMax(Appgen3);
% inja fahmidam ke manfi haye extreme nazdike mosbat haye extreme hast
% va manfi ra mosbat karde va ba bozorgtarin adade mojaver jam kardim

gendet6vert=genstan6vert.*Appgen6;
gendet6hor=genstan6hor.*Appgen6;
gendet6diag=genstan6diag.*Appgen6;
Appgen5=idwt2(Appgen6,gendet6hor,gendet6vert,gendet6diag,'haar');
% Appgen2(isnan(Appgen2))=0;
Appgen5=RemoveNeg(Appgen5);
Appgen5=ArrangeMax(Appgen5);

gendet5diag=genstan5diag.*Appgen5;
gendet5hor=genstan5hor.*Appgen5;
gendet5vert=genstan5vert.*Appgen5;
Appgen4=idwt2(Appgen5,gendet5hor,gendet5vert,gendet5diag,'haar');
% Appgen1(isnan(Appgen1))=0;
Appgen4=RemoveNeg(Appgen4);
Appgen4=ArrangeMax(Appgen4);

gendet4vert=genstan4vert.*Appgen4;
gendet4diag=genstan4diag.*Appgen4;
gendet4hor=genstan4hor.*Appgen4;
Appgen3=idwt2(Appgen4,gendet4hor,gendet4vert,gendet4diag,'haar');
% Sim(isnan(Sim))=0;
Appgen3=RemoveNeg(Appgen3);
Appgen3=ArrangeMax(Appgen3);

gendet3vert=genstan3vert.*Appgen3;
gendet3diag=genstan3diag.*Appgen3;
gendet3hor=genstan3hor.*Appgen3;
Appgen2=idwt2(Appgen3,gendet3hor,gendet3vert,gendet3diag,'haar');
% Sim(isnan(Sim))=0;
Appgen2=RemoveNeg(Appgen2);
Appgen2=ArrangeMax(Appgen2);

gendet2vert=genstan2vert.*Appgen2;
gendet2diag=genstan2diag.*Appgen2;
gendet2hor=genstan2hor.*Appgen2;
Appgen1=idwt2(Appgen2,gendet2hor,gendet2vert,gendet2diag,'haar');
% Sim(isnan(Sim))=0;
Appgen1=RemoveNeg(Appgen1);
Appgen1=ArrangeMax(Appgen1);

gendet1vert=genstan1vert.*Appgen1;
gendet1diag=genstan1diag.*Appgen1;
gendet1hor=genstan1hor.*Appgen1;
Sim=idwt2(Appgen1,gendet1hor,gendet1vert,gendet1diag,'haar');
% Sim(isnan(Sim))=0;
Sim=RemoveNeg(Sim);
Sim=ArrangeMax(Sim);
%% ColorMap

A=zeros(11,1);
B=ones(11,1);
C=linspace(0,1,11)';
D=linspace(1,0,11)';
E=linspace(1,0,6)';
F=ones(6,1);
G=zeros(6,1);
H=linspace(1,0.5,6)';
FirstCol=[E;A;A;C;B;H];
SecondCol=[E;C;B;B;D;G];
ThirdCol=[F;B;D;A;A;G];
ALL=[FirstCol,SecondCol,ThirdCol];
ALL(6,:)=[];
ALL(16,:)=[];
ALL(26,:)=[];
ALL(36,:)=[];
ALL(46,:)=[];

% R_Sim=DC(R,Sim)
% Gen1_Obs1=DC(cA1,Appgen1)
% Gen2_Obs2=DC(cA2,Appgen2)
% Gen3_Obs3=DC(cA3,Appgen3)
% figure; 
% subplot(1,2,1); imagesc(R);
% colorbar
% subplot(1,2,2); imagesc(Sim);
% colorbar
% colormap(ALL)
% figure; 
% subplot(1,2,1); imagesc(cA1);
% colorbar
% subplot(1,2,2); imagesc(Appgen1);
% colorbar
% colormap(ALL)
% figure; 
% subplot(1,2,1); imagesc(cA2);
% colorbar
% subplot(1,2,2); imagesc(Appgen2);
% colorbar
% colormap(ALL)
% figure; 
% subplot(1,2,1); imagesc(cA3);
% colorbar
% subplot(1,2,2); imagesc(Appgen3);
% colorbar
% colormap(ALL)
% figure; 
% subplot(1,2,1); imagesc(cA4);
% colorbar
% subplot(1,2,2); imagesc(Appgen4);
% colorbar
% colormap(ALL)
% figure; 
% subplot(1,2,1); imagesc(cA5);
% colorbar
% subplot(1,2,2); imagesc(Appgen5);
% colorbar
% colormap(ALL)
% figure; 
% subplot(1,2,1); imagesc(cA6);
% colorbar
% subplot(1,2,2); imagesc(Appgen6);
% colorbar
% colormap(ALL)
% newmap = jet; then make the first row 1,1,1 then colormap(newmap);