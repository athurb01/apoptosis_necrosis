function dydt=pysb_D2FC_ODE(t,y,TR);
Comp = 1;
cell_vol= 2.6999999999999998e-12;
nuc_vol= 6.2790697674418602e-13;
cyt_vol= 2.0720930232558138e-12;

k = 24579;
h = 2;

kact_IKK = 0;
kinact_IKK= 0.0030000000000000001;
kneut_IKK= 0.00059999999999999995;
kbind_IkBa_NF= 8.3027000000000003e-19;
kdis_IkBa_NF= 0.050000000000000003;
kphos_IkBa= 1.2288e-19;
kphos_IkBaNF= 6.1439999999999998e-19;
kdeg_pIkBa= 0.10000000000000001;
kimp_NFkB= 0.0025999999999999999;
kexp_NFkB= 0.00017159999999999997;
kexp_IkBaNF= 0.033000000000000002;
kimp_IkBa= 0.00067000000000000002;
kexp_IkBa= 0.0011054999999999999;
ksynth_IkBa= 0.5;
kdeg_tIkBa= 0.00029999999999999997;
kdeg_IkBa= 0.00050000000000000001;
kdeg_bIkBa= 2.1999999999999999e-05;
ksynth_A20= 0.5;
kdeg_tA20= 0.00040000000000000002;
kdeg_A20= 0.0044999999999999997;
kdeg_tTarget= 0.00040000000000000002;
ksynth_Comp= 0.5;
kdeg_tComp= 4.0000000000000003e-05;
kdeg_Comp= 0.00050000000000000001;

NFkB         = y(1);
IkBa         = y(2);
IkBaNFkB     = y(3);
nNFkB        = y(4);
nIkBa        = y(5);
nIkBaNFkB    = y(6);
tIkBa        = y(7);
IKKn         = y(8);
IKK          = y(9);
IKKIkBa      = y(10);
IKKIkBaNFkB  = y(11);
IKKi         = y(12);
tA20         = y(13);
A20          = y(14);
pIkBa        = y(15);
pIkBaNFkB    = y(16);
tInducedTarget = y(17);
tCompetitor  = y(18);
Competitor   = y(19);
dydt         = zeros(19,1);

% Included TRA20 to allow for kbA20 independence of TR
% i.e. TRA20 = TR (is final model used) TRA20=1 gives A20 TNFa independence
TRA20 = TR;

% System ODE's - No IKK complexes %
dydt(1)=kdis_IkBa_NF*IkBaNFkB - kbind_IkBa_NF*IkBa*NFkB - kimp_NFkB*NFkB + kdeg_bIkBa*IkBaNFkB + kexp_NFkB*nNFkB + kdeg_pIkBa*pIkBaNFkB; % Free Cytoplasmic NFkB      
dydt(2)=kdis_IkBa_NF*IkBaNFkB - kbind_IkBa_NF*IkBa*NFkB - kimp_IkBa*IkBa + kexp_IkBa*nIkBa - kdeg_IkBa*IkBa + ksynth_IkBa*tIkBa - kphos_IkBa*IKK*IkBa; % Free Cytoplasmic IkBa               
dydt(3)=kbind_IkBa_NF*IkBa*NFkB - kdis_IkBa_NF*IkBaNFkB + kexp_IkBaNF*nIkBaNFkB - kdeg_bIkBa*IkBaNFkB - kphos_IkBaNF*IKK*IkBaNFkB; % Cytoplasmic NFkB-IkBa            
dydt(4)=kdis_IkBa_NF*nIkBaNFkB - kbind_IkBa_NF*nIkBa*nNFkB + kimp_NFkB*NFkB - kexp_NFkB*nNFkB;             % Free Nuclear NFkB   
dydt(5)=kdis_IkBa_NF*nIkBaNFkB - kbind_IkBa_NF*nIkBa*nNFkB + kimp_IkBa*IkBa - kexp_IkBa*nIkBa - kdeg_IkBa*nIkBa;   % Free Nuclear IkBa
dydt(6)=kbind_IkBa_NF*nIkBa*nNFkB - kdis_IkBa_NF*nIkBaNFkB - kexp_IkBaNF*nIkBaNFkB;          % Nuclear NFkB-IkBa  
dydt(7)=0.174698*((nNFkB/k)^h/((nNFkB/k)^h + 1)) - kdeg_tIkBa*tIkBa;  % tIkBa 
dydt(8)=kneut_IKK*IKKi - TR*kact_IKK*IKKn;  % IKKn
dydt(9)=(TR*kact_IKK)*IKKn - kinact_IKK*IKK;  % IKK - Free active IKK          
dydt(10)=0;   % IKKIkBa - IkB bound active IKK
dydt(11)=0;   % IKKIkBaNFkB  trimeric IkBa-NFKB bound active IKK
dydt(12)=kinact_IKK*IKK - kneut_IKK*IKKi;    % IKKi - Inactive IKK
dydt(13)=0.249569*((nNFkB/k)^(h+1))/((nNFkB/k)^(h+1) + (Competitor/k)^(h+1) + 1) - kdeg_tA20*tA20;      %  tA20 - A20 transript   
dydt(14)=ksynth_A20*tA20 - kdeg_A20*A20;  % A20
dydt(15)=kphos_IkBa*IKK*IkBa - kdeg_pIkBa*pIkBa;  % pIkBa - Phosphorylated IkBa
dydt(16)=kphos_IkBaNF*IKK*IkBaNFkB - kdeg_pIkBa*pIkBaNFkB; % pIkBaNFkB - NFkB bound phosphoryated IkBa
dydt(17)=0.249569*((nNFkB/k)^(h+1))/((nNFkB/k)^(h+1) + (Competitor/(k/2))^(h+1) + 1) - kdeg_tTarget*tInducedTarget;    % tInducedTarget transcription
if Comp == 1
    dydt(18)=0.174698*((nNFkB/k)^(h+1)/((nNFkB/k)^(h+1)  + (Competitor/(k))^(h+1)+ 1)) - kdeg_tComp*tCompetitor;%**tCompetitor competitor 
    dydt(19)=ksynth_Comp*tCompetitor - kdeg_Comp*Competitor;% Competitor protein

end