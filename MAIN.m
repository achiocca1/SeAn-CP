%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directoryRES = 'Results_2D_TE_TO_NON_PROP_2'; % Directory where all the ANSYS filea are located
blocklength = 40; % How many load blocks
times_bb = [36,40]; % Bending Blocks
NameFile = 'Prova.txt';
kFS  = 0.4; % Fatemi Socie
Sy = 300; % Fatemi Socie
INF = 0.01;% Infittimento
symm = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% START DATA READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Files = dir(directoryRES); 
% Read Stress Strain from ANSYS files
FileNames = [];
for K=1:length(Files)
    
    if Files(K).bytes == 0
        continue % Eliminate the files with zero bytes extension
    else
        file = Files(K).name;
        FileNames = [FileNames, string(file)];
    end
    
end
% Import data from files
index = 1; % Necessary for good preallocation of matrix seizes
f = waitbar(0,'Simulation in progress...'); % PROGRESS BAR Option
for i = 1 : length(FileNames)
    
    file = csvread(strcat(directoryRES,"\",FileNames(i)),2,1);
    nodenumber = sscanf(FileNames(i), 'risultati_%d.csv'); % get the node number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% FINISH DATA READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Stress and Strain matrix creation trough cell array
    % Preallocate cell memory
    E = cell(1, blocklength);
    S = cell(1, blocklength);
    EigenvalE = cell(1, blocklength);
    EigenvectE = cell(1, blocklength);
    EigenvalS = cell(1, blocklength);
    EigenvectS = cell(1, blocklength);
    GammaMax = cell(1, blocklength);
    for j = 1 : blocklength
        E0 = [file(j,1) + file(j,7), (file(j,4) + file(j,10))/2, (file(j,6) + file(j,12))/2;
            (file(j,4) + file(j,10))/2, file(j,2) + file(j,8), (file(j,5) + file(j,11))/2;
            (file(j,6) + file(j,12))/2, (file(j,5) + file(j,11))/2, file(j,3) + file(j,9)];
        E{j} = E0; % Strain Tensor
    end
    
    for j = 1 : blocklength
        S0 = [file(j,13), file(j,16), file(j,18);
            file(j,16), file(j,14), file(j,17);
            file(j,18), file(j,17), file(j,15)];
        S{j} = S0; % Strain Tensor
    end
    

end
    index = index + 1;
    
    %% PROGRESS BAR
    waitbar(i/length(FileNames))
    
    
%Delta eps    
DeltaE = E{1, times_bb(1)} - E{1, times_bb(2)};
% Delta Sigma
DeltaS = S{1, times_bb(1)} - S{1, times_bb(2)};

% Calculate and sort principal strains
[V0, D0] = eig(DeltaE, 'vector');
[D0, ind] = sort(D0, 'descend');
V0 = V0(:, ind);
% Calculate and sort principal stresses S1
[W01, G01] = eig(S{1, times_bb(1)}, 'vector');
[G01, ind] = sort(G01, 'descend');
W01 = W01(:, ind);
% Calculate and sort principal stresses S2
[W02, G02] = eig(S{1, times_bb(2)}, 'vector');
[G02, ind] = sort(G02, 'descend');
W02 = W02(:, ind);

EigenvectE = V0; % Strain Tensor
EigenvalE = D0; % Strain Tensor
EigenvectS1 = W01; % Stress Tensor 1
EigenvalS1 = G01; % Stress Tensor 1
EigenvectS2 = W02; % Stress Tensor 2
EigenvalS2 = G02; % Stress Tensor 2


% Calcultate the rotaion angles XYZ fro Deltaeps
THETA_epsPos = [];
PSI_epsPos = [];
THETA_epsNeg = [];
PSI_epsNeg = [];
ALPHA_eps = [];
FS_eps = [];
% Sigman1_eps = V0*S{1, times_bb(1)}*V0';
% Sigman2_eps = V0*S{1, times_bb(2)}*V0';
%% Track point of Mohr's circle %%
Sigma_1_1 = [];
Tau_1_1 = [];
Sigma_1_2 = [];
Tau_1_2 = [];
DELTA_GAMMA_1 = [];
DELTA_EPS_1 = [];
tic
for i = 0 : INF : pi/4
    % i = i
    MatrPos = V0*[cos(i)   0   sin(i);
        0        1        0;
        -sin(i)   0   cos(i)];

    DeltaGamma = ((EigenvalE(1) - EigenvalE(3))/2)*sin(2*i);
    DELTA_GAMMA_1 = [DELTA_GAMMA_1, DeltaGamma];
    DeltaEps = 0.5*(EigenvalE(1) + EigenvalE(3)) + ((EigenvalE(1) - EigenvalE(3))/2)*cos(2*i);
    DELTA_EPS_1 = [DELTA_EPS_1, DeltaEps];
    
    SIGMA = MatrPos'*S{1, times_bb(1)}*MatrPos;
    Sigman11 = SIGMA(1,1);

    SIGMA = MatrPos'*S{1, times_bb(2)}*MatrPos;
    Sigman21 = SIGMA(1,1);
    

    FS1 = DeltaGamma.*(1 + kFS.*(Sigman11./Sy));
    FS2 = DeltaGamma.*(1 + kFS.*(Sigman21./Sy));

    % i = i + pi/2
    MatrPos = V0*[cos(i + pi/2)   0   sin(i + pi/2);
        0        1        0;
        -sin(i + pi/2)   0   cos(i + pi/2)];

    SIGMA = MatrPos'*S{1, times_bb(1)}*MatrPos;
    Sigman11 = SIGMA(1,1);

    SIGMA = MatrPos'*S{1, times_bb(2)}*MatrPos;
    Sigman21 = SIGMA(1,1);
    
    FS3 = DeltaGamma.*(1 + kFS.*(Sigman11./Sy));
    FS4 = DeltaGamma.*(1 + kFS.*(Sigman21./Sy));

    FS_eps = [FS_eps, max([FS1, FS2, FS3, FS4])];
    
    if max(FS1, FS2) > max(FS2, FS3)
        q = i + pi/2 + symm;
    else
        q = i + symm;
    end

        MatrPos = V0*[cos(q)   0   sin(q);
            0        1        0;
            -sin(q)   0   cos(q)];
        MatrNeg = V0*[cos(-q)   0   sin(-q);
            0        1        0;
            -sin(-q)   0   cos(-q)];

    ThetaPos = atan2(sqrt(MatrPos(1,3)^2+MatrPos(2,3)^2), MatrPos(3,3));
    PsiPos = atan2(MatrPos(2,3), MatrPos(1,3));
    ThetaNeg = atan2(sqrt(MatrNeg(1,3)^2+MatrNeg(2,3)^2), MatrNeg(3,3));
    PsiNeg = atan2(MatrNeg(2,3), MatrNeg(1,3));
    THETA_epsPos = [THETA_epsPos, ThetaPos];
    PSI_epsPos = [PSI_epsPos, PsiPos];
    THETA_epsNeg = [THETA_epsNeg, ThetaNeg];
    PSI_epsNeg = [PSI_epsNeg, PsiNeg];
    %% Track point of Mohr's circle %%  
    % First normal direction z
    n = MatrPos(:,3); % Normal
    Sig = S{1, times_bb(1)}; % Stress tensor
    Sigma_n = Sig(1,1)*n(1)^2 + Sig(2,2)*n(2)^2 + Sig(3,3)*n(3)^2 + 2*(Sig(1,2)*n(1)*n(2) + Sig(2,3)*n(2)*n(3) + Sig(1,3)*n(1)*n(3));
    T_1 = Sig(1,1)*n(1) + Sig(1,2)*n(2) + Sig(1,3)*n(3);
    T_2 = Sig(2,1)*n(1) + Sig(2,2)*n(2) + Sig(2,3)*n(3);
    T_3 = Sig(3,1)*n(1) + Sig(3,2)*n(2) + Sig(3,3)*n(3);
    T_squar = T_1^2 + T_2^2 + T_3^2;
    Tau_n = (T_squar - Sigma_n^2)^0.5;
    
    Sigma_1_1 = [Sigma_1_1, Sigma_n];
    Tau_1_1 = [Tau_1_1, Tau_n];
    % Second normal direction x
    n = MatrPos(:,3); % Normal
    Sig = S{1, times_bb(2)}; % Stress tensor
    Sigma_n = Sig(1,1)*n(1)^2 + Sig(2,2)*n(2)^2 + Sig(3,3)*n(3)^2 + 2*(Sig(1,2)*n(1)*n(2) + Sig(2,3)*n(2)*n(3) + Sig(1,3)*n(1)*n(3));
    T_1 = Sig(1,1)*n(1) + Sig(1,2)*n(2) + Sig(1,3)*n(3);
    T_2 = Sig(2,1)*n(1) + Sig(2,2)*n(2) + Sig(2,3)*n(3);
    T_3 = Sig(3,1)*n(1) + Sig(3,2)*n(2) + Sig(3,3)*n(3);
    T_squar = T_1^2 + T_2^2 + T_3^2;
    Tau_n = (T_squar - Sigma_n^2)^0.5;
    
    Sigma_1_2 = [Sigma_1_2, Sigma_n];
    Tau_1_2 = [Tau_1_2, Tau_n];
    % Track point of Mohr's circle %
end

% Calcultate the rotaion angles XYZ for S1
THETA_S1Pos = [];
PSI_S1Pos = [];
THETA_S1Neg = [];
PSI_S1Neg = [];
ALPHA_S1 = [];
FS_S1 = [];
%DeltaE_S1 = W01*DeltaE*W01';
%% Track point of Mohr's circle %%
Sigma_2_1 = [];
Tau_2_1 = [];
Sigma_2_2 = [];
Tau_2_2 = [];
DELTA_GAMMA_2 = [];
DELTA_EPS_2 = [];
for i = 0 : INF : pi/4
    MatrPos = W01*[cos(i + pi/2 + symm)   0   sin(i + pi/2 + symm);
        0        1        0;
        -sin(i + pi/2 + symm)   0   cos(i + pi/2 + symm)];
    MatrNeg = W01*[cos(-(i + pi/2 + symm))   0   sin(-(i + pi/2 + symm));
        0        1        0;
        -sin(-(i + pi/2 + symm))   0   cos(-(i + pi/2 + symm))];

    ThetaPos = atan2(sqrt(MatrPos(1,3)^2+MatrPos(2,3)^2), MatrPos(3,3));
    PsiPos = atan2(MatrPos(2,3), MatrPos(1,3));
    ThetaNeg = atan2(sqrt(MatrNeg(1,3)^2+MatrNeg(2,3)^2), MatrNeg(3,3));
    PsiNeg = atan2(MatrNeg(2,3), MatrNeg(1,3));
    
    THETA_S1Pos = [THETA_S1Pos, ThetaPos];
    PSI_S1Pos = [PSI_S1Pos, PsiPos];
    THETA_S1Neg = [THETA_S1Neg, ThetaNeg];
    PSI_S1Neg = [PSI_S1Neg, PsiNeg];


    Sigman11 = ((EigenvalS1(1) + EigenvalS1(3))/2) + ((EigenvalS1(1) - EigenvalS1(3))/2)*cos(2*i);
    Sigma_2_1 = [Sigma_2_1, Sigman11];
    Tau11 = 0.5*(EigenvalS1(1) - EigenvalS1(3))*sin(2*i);
    Tau_2_1 = [Tau_2_1, Tau11];
    
    DeltaE_S1_1 = MatrPos'*DeltaE*MatrPos;
    DeltaGamma = norm(DeltaE_S1_1(2:3,1));
    
    FS1 = DeltaGamma.*(1 + kFS.*(Sigman11./Sy));

    SIGMA = MatrPos'*S{1, times_bb(2)}*MatrPos;
    Sigman21 = SIGMA(3,3);

    FS2 = DeltaGamma.*(1 + kFS.*(Sigman21./Sy));

    FS_S1 = [FS_S1, max([FS1, FS2])];
    %% Track point of Mohr's circle %%  
    % First normal direction z
    n = MatrPos(:,3); % Normal
    Sig = DeltaE; % Strain range tensor
    Sigma_n = Sig(1,1)*n(1)^2 + Sig(2,2)*n(2)^2 + Sig(3,3)*n(3)^2 + 2*(Sig(1,2)*n(1)*n(2) + Sig(2,3)*n(2)*n(3) + Sig(1,3)*n(1)*n(3));
    T_1 = Sig(1,1)*n(1) + Sig(1,2)*n(2) + Sig(1,3)*n(3);
    T_2 = Sig(2,1)*n(1) + Sig(2,2)*n(2) + Sig(2,3)*n(3);
    T_3 = Sig(3,1)*n(1) + Sig(3,2)*n(2) + Sig(3,3)*n(3);
    T_squar = T_1^2 + T_2^2 + T_3^2;
    Tau_n = (T_squar - Sigma_n^2)^0.5;
    
    DELTA_EPS_2 = [DELTA_EPS_2, Sigma_n];
    DELTA_GAMMA_2 = [DELTA_GAMMA_2, Tau_n];
    % Second normal direction x
    n = MatrPos(:,3); % Normal
    Sig = S{1, times_bb(2)}; % Stress tensor
    Sigma_n = Sig(1,1)*n(1)^2 + Sig(2,2)*n(2)^2 + Sig(3,3)*n(3)^2 + 2*(Sig(1,2)*n(1)*n(2) + Sig(2,3)*n(2)*n(3) + Sig(1,3)*n(1)*n(3));
    T_1 = Sig(1,1)*n(1) + Sig(1,2)*n(2) + Sig(1,3)*n(3);
    T_2 = Sig(2,1)*n(1) + Sig(2,2)*n(2) + Sig(2,3)*n(3);
    T_3 = Sig(3,1)*n(1) + Sig(3,2)*n(2) + Sig(3,3)*n(3);
    T_squar = T_1^2 + T_2^2 + T_3^2;
    Tau_n = (T_squar - Sigma_n^2)^0.5;
    
    Sigma_2_2 = [Sigma_2_2, Sigma_n];
    Tau_2_2 = [Tau_2_2, Tau_n];
    % Track point of Mohr's circle %
end

% Calcultate the rotaion angles XYZ for S2
THETA_S2Pos = [];
PSI_S2Pos = [];
THETA_S2Neg = [];
PSI_S2Neg = [];
ALPHA_S2 = [];
FS_S2 = [];
%DeltaE_S2 = W02*DeltaE*W02';
%% Track point of Mohr's circle %%
Sigma_3_1 = [];
Tau_3_1 = [];
Sigma_3_2 = [];
Tau_3_2 = [];
DELTA_GAMMA_3 = [];
DELTA_EPS_3 = [];
for i = 0 : INF : pi/4
    MatrPos = W02*[cos(i + pi/2 + symm)   0   sin(i + pi/2 + symm);
        0        1        0;
        -sin(i + pi/2 + symm)   0   cos(i + pi/2 + symm)];
    MatrNeg = W02*[cos(-(i + pi/2 + symm))   0   sin(-(i + pi/2 + symm));
        0        1        0;
        -sin(-(i + pi/2 + symm))   0   cos(-(i + pi/2 + symm))];

    ThetaPos = atan2(sqrt(MatrPos(1,3)^2+MatrPos(2,3)^2), MatrPos(3,3));
    PsiPos = atan2(MatrPos(2,3), MatrPos(1,3));
    ThetaNeg = atan2(sqrt(MatrNeg(1,3)^2+MatrNeg(2,3)^2), MatrNeg(3,3));
    PsiNeg = atan2(MatrNeg(2,3), MatrNeg(1,3));
    
    THETA_S2Pos = [THETA_S2Pos, ThetaPos];
    PSI_S2Pos = [PSI_S2Pos, PsiPos];
    THETA_S2Neg = [THETA_S2Neg, ThetaNeg];
    PSI_S2Neg = [PSI_S2Neg, PsiNeg];

    Sigman21 = ((EigenvalS2(1) + EigenvalS2(3))/2) + ((EigenvalS2(1) - EigenvalS2(3))/2)*cos(2*i);
    Sigma_3_2 = [Sigma_3_2, Sigman21];
    Tau21 = 0.5*(EigenvalS2(1) - EigenvalS2(3))*sin(2*i);
    Tau_3_2 = [Tau_3_2, Tau21];
    
    DeltaE_S2_1 = MatrPos'*DeltaE*MatrPos;
    DGamma = 0;
    DeltaGamma = norm(DeltaE_S2_1(2:3,1));
    
    FS2 = DeltaGamma.*(1 + kFS.*(Sigman21./Sy));
   
    FS_S2 = [FS_S2, FS2];
    
    %% Track point of Mohr's circle %%   
    % First normal direction z
    n = MatrPos(:,3); % Normal
    Sig = DeltaE; % Strain range tensor
    Sigma_n = Sig(1,1)*n(1)^2 + Sig(2,2)*n(2)^2 + Sig(3,3)*n(3)^2 + 2*(Sig(1,2)*n(1)*n(2) + Sig(2,3)*n(2)*n(3) + Sig(1,3)*n(1)*n(3));
    T_1 = Sig(1,1)*n(1) + Sig(1,2)*n(2) + Sig(1,3)*n(3);
    T_2 = Sig(2,1)*n(1) + Sig(2,2)*n(2) + Sig(2,3)*n(3);
    T_3 = Sig(3,1)*n(1) + Sig(3,2)*n(2) + Sig(3,3)*n(3);
    T_squar = T_1^2 + T_2^2 + T_3^2;
    Tau_n = (T_squar - Sigma_n^2)^0.5;
    
    DELTA_EPS_3 = [DELTA_EPS_3, Sigma_n];
    DELTA_GAMMA_3 = [DELTA_GAMMA_3, Tau_n];
    % Second normal direction x
    n = MatrPos(:,3); % Normal
    Sig = S{1, times_bb(1)}; % Stress tensor
    Sigma_n = Sig(1,1)*n(1)^2 + Sig(2,2)*n(2)^2 + Sig(3,3)*n(3)^2 + 2*(Sig(1,2)*n(1)*n(2) + Sig(2,3)*n(2)*n(3) + Sig(1,3)*n(1)*n(3));
    T_1 = Sig(1,1)*n(1) + Sig(1,2)*n(2) + Sig(1,3)*n(3);
    T_2 = Sig(2,1)*n(1) + Sig(2,2)*n(2) + Sig(2,3)*n(3);
    T_3 = Sig(3,1)*n(1) + Sig(3,2)*n(2) + Sig(3,3)*n(3);
    T_squar = T_1^2 + T_2^2 + T_3^2;
    Tau_n = (T_squar - Sigma_n^2)^0.5;
    
    Sigma_3_1 = [Sigma_3_1, Sigma_n];
    Tau_3_1 = [Tau_3_1, Tau_n];
    % Track point of Mohr's circle %
end
toc

THETAPos = [THETA_epsPos, THETA_S1Pos, THETA_S2Pos];
PSIPos = [PSI_epsPos, PSI_S1Pos, PSI_S2Pos];
THETANeg = [THETA_epsNeg, THETA_S1Neg, THETA_S2Neg];
PSINeg = [PSI_epsNeg, PSI_S1Neg, PSI_S2Neg];
[FSmax, FSi] = max([FS_eps, FS_S1, FS_S2]);
VectSIGMA_n_EPS = [DELTA_EPS_1, DELTA_EPS_2, DELTA_EPS_3];
VectSIGMA_n_S1 = [Sigma_1_1, Sigma_2_1, Sigma_3_1];
VectSIGMA_n_S2 = [Sigma_1_2, Sigma_2_2, Sigma_3_2];
VectTAU_n_EPS = [DELTA_GAMMA_1, DELTA_GAMMA_2, DELTA_GAMMA_3];
VectTAU_n_S1 = [Tau_1_1, Tau_2_1, Tau_3_1];
VectTAU_n_S2 = [Tau_1_2, Tau_2_2, Tau_3_2];
FSmax
THETAPos(FSi)
PSIPos(FSi)
THETANeg(FSi)
PSINeg(FSi)
delete(f) % wait bar deleted

save Var.mat THETA_epsPos PSI_epsPos FS_eps THETA_S1Pos PSI_S1Pos FS_S1 THETA_S2Pos PSI_S2Pos FS_S2 THETA_epsNeg PSI_epsNeg THETA_S1Neg PSI_S1Neg THETA_S2Neg PSI_S2Neg