function spectrumSD = semispectrum_2DSD_GenPatt(sdisc,doplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script computes the eiegevalues of the Spectral Difference 
%% operator for the 2 D linear advection equation using the wave 
%% propagation algorithm. One generating pattern is used. 
%%
%% In the algorithm no particular position of the solution points is 
%% assumed. This additional "degree of freedom" allows to study the 
%% effect of the position of the solution points on the time marching 
%% scheme. 
%%
%%
%%
%% Author:              Matteo Parsani
%% Date:                2011-05-24
%% Last Modified Date:  2012-01-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dereference some parameters 
order = sdisc.order;
upwindPar = sdisc.upwindPar;
KStep = sdisc.KStep;
thetaStep = sdisc.thetaStep;
psiStep = sdisc.psiStep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D Flux and solution points position. 
% The coordinates of the flux points are based on the stability accuracy 
% properties of the scheme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch order
    case 1
        % 1st order
        fluxPoints1D = [ -1.0 ; 1.0 ];
        solPoints1D = [ 0 ];
    case 2
        % 2nd order
        fluxPoints1D = [ -1.0 ; 0.0 ; 1.0 ];
        solPoints1D = [ -1.0 ;       1.0 ];
    case 3
        % 3rd order
        fluxPoints1D = [ -1.0 ; -0.58 ; 0.58 ; 1.0 ];
        solPoints1D = [ -1.0 ;         0.58 ; 1.0 ];
    case 4
        % 4th order
        fluxPoints1D = [ -1.0 ; -0.78 ; 0.0 ; 0.78 ; 1.0 ];
        solPoints1D = [ -1.0 ; -0.78 ;       0.78 ; 1.0 ];
    case 5
        % 5th order
        fluxPoints1D = [ -1.0 ; -0.83 ; -0.36 ; 0.36 ; 0.83 ; 1.0 ];
        solPoints1D = [ -1.0 ; -0.83 ;         0.36 ; 0.83 ; 1.0 ];
    case 6
        % 6th order
        fluxPoints1D = [ -1.0 ; -0.88 ; -0.53 ; 0.0 ; 0.53 ; 0.88 ; 1.0 ];
        solPoints1D = [ -1.0 ; -0.88 ; -0.53 ;       0.53 ; 0.88 ; 1.0 ];
    otherwise
        disp('Error: min order 1, max order 6');
end


nbrFluxPoints1D = length(fluxPoints1D);
nbrSolPoints1D = length(solPoints1D);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct 1D matrices.
% This matrices are used later as a building blocks for the 2D matrices.
% They are independent of the direction of the convective velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix for the calculation of the solution to the flux points
extrSolToFlux1D = ones(nbrFluxPoints1D, nbrSolPoints1D);

for iFlux = 1:nbrFluxPoints1D
    for iSol = 1:nbrSolPoints1D
        extrSolToFlux1D(iFlux,iSol) = 1.0;
        for iFactor = 1:nbrSolPoints1D
            if iFactor ~= iSol
                extrSolToFlux1D(iFlux,iSol) = extrSolToFlux1D(iFlux,iSol)*((fluxPoints1D(iFlux)-solPoints1D(iFactor))/(solPoints1D(iSol)-solPoints1D(iFactor)));
            end
        end
    end
end

extrSolToFlux1D;


% Matrix for the calculation of the flux to the solution points
extrFluxToSol1D = ones(nbrSolPoints1D, nbrFluxPoints1D);

for iSol = 1:nbrSolPoints1D
    for iFlux = 1:nbrFluxPoints1D
        for iFactor = 1:nbrFluxPoints1D
            if iFactor ~= iFlux
                extrFluxToSol1D(iSol,iFlux) = extrFluxToSol1D(iSol,iFlux)*((solPoints1D(iSol)-fluxPoints1D(iFactor))/(fluxPoints1D(iFlux)-fluxPoints1D(iFactor)));
            end
        end
    end
end

extrFluxToSol1D;


% Matrix for the calculation of the flux derivative at the solution points
derivFluxInSolPoints1D = zeros(nbrSolPoints1D,nbrFluxPoints1D);
term = 0;
for iFlux = 1:nbrFluxPoints1D 
    for iSol = 1:nbrSolPoints1D
        for iTerm = 1:nbrFluxPoints1D
            if iTerm ~= iFlux
                term = 1.0/(fluxPoints1D(iFlux)-fluxPoints1D(iTerm));
                for iFactor = 1:nbrFluxPoints1D
                    if (iFactor ~= iTerm) && (iFactor ~= iFlux)
                        term = term*(solPoints1D(iSol)-fluxPoints1D(iFactor))/(fluxPoints1D(iFlux)-fluxPoints1D(iFactor));
                    end
                end
                derivFluxInSolPoints1D(iSol,iFlux) = derivFluxInSolPoints1D(iSol,iFlux) + term;
            end
        end
    end
end

derivFluxInSolPoints1D;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating pattern of the quadrilateral (2D) spectral difference cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating pattern: square cell 
% Nodes numbered in counter-clock-wise
% First point : (0.0,0.0)
% Second point: (1.0,0.0) 
% Third point : (1.0,1.0) 
% Fourth point: (0.0,1.0)

genPatNodes = [0, 0; 1, 0; 1, 1; 0, 1];


% Face size
faceSize = genPatNodes(2,1) - genPatNodes(1,1);

% Transformation matrix, i.e. {x_vec} = [T] {csi_vec}
T = zeros(2,2);
T =[
    (genPatNodes(2, 1) - genPatNodes(1, 1)), (genPatNodes(4, 1) - genPatNodes(1, 1)); 
    (genPatNodes(2, 2) - genPatNodes(1, 2)), (genPatNodes(4, 2) - genPatNodes(1, 2));     
   ];

T;

% Inverse of the transformation matrix
invT = inv(T);


% Normals to generating pattern faces
% 1st normal is the normal to the vertical faces (to compute X-fluxes)
% 2nd normal is the normal to the horizontal faces (to compute Y-fluxes)
genPatFaceNormals = zeros(2,2);
genPatFaceNormals(1,1) = +(genPatNodes(4, 2) - genPatNodes(1, 2))/faceSize;
genPatFaceNormals(1,2) = -(genPatNodes(4, 1) - genPatNodes(1, 1))/faceSize;
genPatFaceNormals(2,1) = -(genPatNodes(2, 2) - genPatNodes(1, 2))/faceSize;
genPatFaceNormals(2,2) = +(genPatNodes(2, 1) - genPatNodes(1, 1))/faceSize;

genPatFaceNormals;

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct 2D matrices.
% This matrices are constructed by using the 1D matrices assembled before.
% They are independent of the direction of the convective velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Direction of the convective speed
%psiStep = 0.1;
psiValues = [0:psiStep:2*pi];
nbrPsi = length(psiValues);


% Define wave number vector. This leads to two for loops, i.e. one for
% the module of the wave number vector (KValues) and one for the 
% orientation (thetaValues).
    
% Wave number module
%KStep = 0.1;
KValues = [0:KStep:2*pi];
nbrK = length(KValues);
    
% Wave number orientation vector
%thetaStep = 0.1;
thetaValues = [0:thetaStep:2*pi];
nbrTheta = length(thetaValues);

% Number of DOF in 2D
nbrSolPoints2D = nbrSolPoints1D*nbrSolPoints1D;

% Array that will contain the spectrum of the 2D linear operator L.
spectrumSD = [];

for iPsi = 1:nbrPsi
    
    psi = psiValues(iPsi);
    cPsi = cos(psi);  
    sPsi = sin(psi);
    
    
    % Matricies for calculation of the solution in the flux points in the 
    % csi direction, i.e. M0, MM1, MP1.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M0 matrix
    M0 = zeros(nbrFluxPoints1D*nbrFluxPoints1D,nbrSolPoints2D);
    for iYFlux = 1: nbrFluxPoints1D 
        for iXFlux = 2:nbrFluxPoints1D-1 
            for iYSol = 1:nbrSolPoints1D
                for iXSol = 1:nbrSolPoints1D
                    M0(nbrFluxPoints1D*(iYFlux-1) + iXFlux, nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                        extrSolToFlux1D(iXFlux, iXSol)*extrSolToFlux1D(iYFlux, iYSol)*(invT(1, 1)*cPsi + invT(1, 2)*sPsi);
                end 
            end
        end
    end
    
    M0;
    
    
    % Riemann flux for the vertical faces
    rSolvFactor = invT(1, 1)*cPsi*genPatFaceNormals(1, 1) + invT(1, 2)*sPsi*genPatFaceNormals(1, 2);
    Ap = (rSolvFactor + upwindPar*abs(rSolvFactor))/2;
    Am = (rSolvFactor - upwindPar*abs(rSolvFactor))/2;
    
    
    % MM1 and M0 matrices (left face flux points)
    MM1 = zeros(nbrFluxPoints1D*nbrFluxPoints1D,nbrSolPoints2D);
    for iYFlux = 1:nbrFluxPoints1D
        for iYSol = 1:nbrSolPoints1D
            for iXSol = 1:nbrSolPoints1D
                M0(nbrFluxPoints1D*(iYFlux-1) + 1,nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                    extrSolToFlux1D(iYFlux,iYSol)*extrSolToFlux1D(1,iXSol)*Am;
                
                MM1(nbrFluxPoints1D*(iYFlux-1) + 1,nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                    extrSolToFlux1D(iYFlux,iYSol)*extrSolToFlux1D(nbrFluxPoints1D,iXSol)*Ap;
            end
        end
    end
    
    M0;
    MM1;
    
    
    % MP1 and M0 matrices (right face fkux points)
    MP1 = zeros(nbrFluxPoints1D*nbrFluxPoints1D,nbrSolPoints2D);
    for iYFlux = 1:nbrFluxPoints1D
        for iYSol = 1:nbrSolPoints1D
            for iXSol = 1:nbrSolPoints1D
                M0(nbrFluxPoints1D*(iYFlux-1) + nbrFluxPoints1D,nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                    extrSolToFlux1D(iYFlux,iYSol)*extrSolToFlux1D(nbrFluxPoints1D,iXSol)*Ap;
                
                MP1(nbrFluxPoints1D*(iYFlux-1) + nbrFluxPoints1D,nbrSolPoints1D*(iYSol - 1) + iXSol) = ...
                    extrSolToFlux1D(iYFlux, iYSol)*extrSolToFlux1D(1, iXSol)*Am;
            end
        end
    end
    
    M0;
    MM1;
    MP1; 
    
                    
    % Derivation matrix for the csi-fluxes
    DCsi = zeros(nbrSolPoints2D,nbrFluxPoints1D*nbrFluxPoints1D);

    for iYFlux = 1:nbrFluxPoints1D
        for iXFlux = 1:nbrFluxPoints1D
            for iYSol = 1:nbrSolPoints1D
                for iXSol = 1:nbrSolPoints1D
                    DCsi(nbrSolPoints1D*(iYSol-1) + iXSol,nbrFluxPoints1D*(iYFlux-1) + iXFlux) = ...
                        extrFluxToSol1D(iYSol, iYFlux)*derivFluxInSolPoints1D(iXSol, iXFlux);
                end
            end
        end
    end
    
    DCsi;
    
   
    % Matricies for calculation of the solution in the flux points in the
    % eta direction, i.e. N0, NM1, NP1.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % N0 matrix
    N0 = zeros(nbrFluxPoints1D*nbrFluxPoints1D,nbrSolPoints2D);
        
    for iXFlux = 1:nbrFluxPoints1D
        for iYFlux = 2:nbrFluxPoints1D-1
            for iXSol = 1:nbrSolPoints1D
                for iYSol = 1:nbrSolPoints1D
                    N0(nbrFluxPoints1D*(iYFlux-1) + iXFlux, nbrSolPoints1D*(iYSol - 1) + iXSol) = ...
                        extrSolToFlux1D(iXFlux, iXSol)*extrSolToFlux1D(iYFlux, iYSol)*(invT(2, 1)*cPsi + invT(2, 2)*sPsi);
                end
            end
        end
    end
    
    N0;
    
    
    % Riemann flux for the horizontal faces
    rSolvFactor = invT(2, 1)*cPsi*genPatFaceNormals(2, 1) + invT(2, 2)*sPsi*genPatFaceNormals(2, 2);
    Ap = (rSolvFactor + upwindPar*abs(rSolvFactor))/2;
    Am = (rSolvFactor - upwindPar*abs(rSolvFactor))/2;
    
      
    % NM1 and N0 matrices (lower face flux points)
    NM1 = zeros(nbrFluxPoints1D*nbrFluxPoints1D,nbrSolPoints2D);
    for iXFlux = 1:nbrFluxPoints1D
        for iXSol = 1:nbrSolPoints1D
            for iYSol = 1:nbrSolPoints1D
                N0(nbrFluxPoints1D*(1-1) + iXFlux, nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                    extrSolToFlux1D(iXFlux, iXSol)*extrSolToFlux1D(1, iYSol)*Am;
                
                NM1(nbrFluxPoints1D*(1-1) + iXFlux, nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                    extrSolToFlux1D(iXFlux, iXSol)*extrSolToFlux1D(nbrFluxPoints1D, iYSol)*Ap;
            end
        end
    end
    
    N0;
    NM1;
    
    
    % NP1 and N0 matrices (upper face flux points)
    NP1 = zeros(nbrFluxPoints1D*nbrFluxPoints1D,nbrSolPoints2D);
    for iXFlux = 1:nbrFluxPoints1D
        for iXSol = 1:nbrSolPoints1D
            for iYSol = 1:nbrSolPoints1D
                N0(nbrFluxPoints1D*(nbrFluxPoints1D-1) + iXFlux, nbrSolPoints1D*(iYSol-1) + iXSol) = ...
                    extrSolToFlux1D(iXFlux, iXSol)*extrSolToFlux1D(nbrFluxPoints1D, iYSol)*Ap;
                
                NP1(nbrFluxPoints1D*(nbrFluxPoints1D-1) + iXFlux, nbrSolPoints1D*(iYSol-1) + iXSol) = ... 
                    extrSolToFlux1D(iXFlux, iXSol)*extrSolToFlux1D(1, iYSol)*Am;
            end
        end
    end
    
    N0;
    NM1;
    NP1;
    
    
    % Derivation matrix for the eta-fluxes
    DEta = zeros(nbrSolPoints2D,nbrFluxPoints1D*nbrFluxPoints1D);
    
    for iYFlux = 1:nbrFluxPoints1D
        for iXFlux = 1:nbrFluxPoints1D
            for iYSol = 1:nbrSolPoints1D
                for iXSol = 1:nbrSolPoints1D
                    DEta(nbrSolPoints1D*(iYSol - 1) + iXSol, nbrFluxPoints1D*(iYFlux - 1) + iXFlux) = ...
                        extrFluxToSol1D(iXSol, iXFlux)*derivFluxInSolPoints1D(iYSol, iYFlux);
                end
            end
        end
    end
    
    DEta;
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct spatial operator L, i.e dU/dt = L*U.
    % The construction is based on the classic Von Neumann analysis.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % General definition
    %                     _     _                                        _
    % dU_(i,j)     |a|   |     |                                          |
    % -------- +  ------*|DCsi*|MM1*U_(i-1,j) + M0*U_(i,j) + MP1*U_(i+1,j)|
    %   dt        DeltaB |_    |_                                        _|
    %            
    %                     _                                        _  _
    %                    |                                          |  |
    %             + DEta*|NM1*U_(i,j-1) + N0*U_(i,j) + NP1*U_(i,j+1)|  | = 0
    %                    |_                                        _| _|
    %
    % where DeltaB is the module of the vector B_(1,vect), which 
    % together with B_(2,vect) form a pair of vectors that completely 
    % defined the generating pattern, i.e. the quadrilaterla cell.
    % |a| is the module of the convective velocity.
    %
    % Here DeltaB      = 1 
    %      |a| = a_mod = 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    % The Von Neumann analysis uses: 
    %            
    %  **   U_(i,j) = U_tilde *exp(I*k_vect*(i*B_(1,vect) + j*B_(2,vect)))
    %
    % where:
    % I: imaginary unit
    %
    % k_vect = |k|*{cos(theta), sin(theta)}^T: "wave number vector" which 
    % defines the direction of the initial plane wave
    %
    % B_(1,vect), B_(2,vect): pair of vector that completely defined the 
    % generating pattern, i.e. the quadrilaterla cell
    % 
    % Therefore, substituting the spatial Fourier wawe of the form ** one
    % gets:
    %
    %                             _     _                                        
    %   d                      
    % ----- U_tilde +  
    %  dt              
    %
    %                  _     _
    %        |a|      |     |                                  
    %  +   -------- * |DCsi*|MM1*exp(-I*K_vect*B_(1,vect,ad) +
    %       DeltaB    |_    |_             
    %                                   
    %                        M0                              + 
    %                                                     _
    %                                                      |
    %                        MP1*exp(I*K_vect*B_(1,vect,ad)| +
    %                                                     _|
    %                  _     _                                        
    %                 |     |                                          
    %          +      |DEta*|NM1*exp(-I*K_vect*B_(2,vect,ad) + 
    %                 |_    |_                                
    %
    %                   
    %                        N0                              +     
    %                                                     _  _
    %                                                      |  |
    %                        NP1*exp(I*K_vect*B_(2,vect,ad)|  | * U_tilde
    %                                                     _| _|
    %                
    % where K_vect = k_vect* DeltaB = K*{cos(theta), sin(theta)}^T.
    %
    %
    % Consequently the linear operator L is defined as:
    %                    
    %                            _     _
    %                  |a|      |     |                                  
    % L(DeltaB, p) = -------- * |DCsi*|MM1*exp(-I*K_vect*B_(1,vect,ad) +
    %                 DeltaB    |_    |_             
    %                                   
    %                                  M0                              + 
    %                                                     _
    %                                                      |
    %                                  MP1*exp(I*K_vect*B_(1,vect,ad)| +
    %                                                     _|
    %                            _     _                                        
    %                           |     |                                          
    %          +                |DEta*|NM1*exp(-I*K_vect*B_(2,vect,ad) + 
    %                           |_    |_                                
    %
    %                   
    %                                  N0                              +     
    %                                                               _  _
    %                                                                |  |
    %                                  NP1*exp(I*K_vect*B_(2,vect,ad)|  |
    %                                                               _| _|
    % 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % Module of the convective velocity
    a_mod = 1;
    
    % Appply the transformation from physical system to reference system 
    % for two points of the physical generating pattern. 
    p1Map= invT*genPatNodes(1,:)';
    
    p2Map = invT*genPatNodes(2,:)';
    
    % Length of the path connecting them, i.e. DeltaB
    DeltaB = sqrt((p1Map(1,1)-p2Map(1,1))^2 + (p1Map(2,1)-p2Map(2,1))^2);
    
    % Component of the nondimensional generating pattern's vectors
    B1_ad(1) = T(1,1);
    B1_ad(2) = T(2,1);
    
    B2_ad(1) = T(1,2);
    B2_ad(2) = T(2,2);
    
    
    % Loop over the wave number module
    for iK = 1:nbrK
        K = KValues(iK);
        
        % Loop over the orientation of the wave number vector
        for iTheta = 1:nbrTheta
            theta = thetaValues(iTheta);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % Nmerical frequencies %
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute matrix for the calculation of the numerical
            % frequencies
            EigSystMatrix = -sqrt(-1)*(2*DCsi*(M0+ MM1*exp(-sqrt(-1)*K*(T(1,1)*cPsi + T(2,1)*sPsi)) + MP1*exp(sqrt(-1)*K*(T(1,1)*cPsi + T(2,1)*sPsi))) + ...
                2*DEta*(N0 + NM1*exp(-sqrt(-1)*K*(T(1,2)*cPsi + T(2,2)*sPsi)) + NP1*exp(sqrt(-1)*K*(T(1,2)*cPsi + T(2,2)*sPsi))));
            
            % Compute eigenvalues == numerical frequencies
            numFreq = eig(EigSystMatrix);
     
            reNumFreq = real(numFreq);
            imagNumFreq = imag(numFreq);
     
            % Check that the real part of the numerical frequencies is
            % smaller than zero. 
            % Note that numerical frequency = (-I*numFreq) 
            for i = 1: length(numFreq)
                if imagNumFreq(i)>0
                    imagNumFreq(i);
                end
            end
    
            % Plot the Fourier footprint
            %hold on;
            %plot(imagNumFreq,reNumFreq,'o');
            
            % Construct operator the 2D linear operator L
            L = -2*a_mod/DeltaB*(DCsi*(MM1*exp(-sqrt(-1)*K*(cos(theta)*B1_ad(1)+sin(theta)*B1_ad(2))) + ...
                                    M0 + ...
                                    MP1*exp(sqrt(-1)*K*(cos(theta)*B1_ad(1)+sin(theta)*B1_ad(2)))) + ...
                              DEta*(NM1*exp(-sqrt(-1)*K*(cos(theta)*B2_ad(1)+sin(theta)*B2_ad(2))) + ...
                                    N0 + ...
                                    NP1*exp(sqrt(-1)*K*(cos(theta)*B2_ad(1)+sin(theta)*B2_ad(2)))));  
            
            % Compute the eigenvalues of the 2D linear operator L
            lambda = eig(L);
            
            % Append new computed eiegenvalues to SpectrumSD.
            % SpectrumSD will be returned to convex_form for the
            % optimization of the RK scheme.
            for iEigen = 1: nbrSolPoints2D
                spectrumSD(length(spectrumSD(:)) + iEigen) = lambda(iEigen);
            end
             
            % Check that the real part of the eigenvalues is smaller than 
            % zero. 
            for i = 1: length(lambda)
                if real(lambda(i))>0
                    lambda(i);
                end
            end
            
            

        end % Close loop of the wave number orientation
        
    end % Close loop of the wave number module
    
end % Close loop of the convective velocity orientation

spectrumSD = spectrumSD';


% Plot
%*****
if doplot
    % Plot eigenvalues for checking their position in the complex plane
    %******************************************************************
    figure;
    plot(real(spectrumSD(:)),imag(spectrumSD(:)),'o','markersize',8, 'MarkerEdgeColor','k',...
                                                                 'MarkerFaceColor','k');

    hand = gca;
    title_handle = title(['SD spectrum for the linear advection equation, order=',int2str(order)]);
    set(title_handle,'FontSize',16);
    set(title_handle,'FontWeight','bold');

    xlabel('Re(eigenvalue)');
    xlab = get(hand, 'xlabel');
    set(xlab,'FontSize',16);
    set(xlab,'FontWeight','bold');

    ylabel('Im(eigenvalue)');
    ylab = get(hand, 'ylabel');
    set(ylab,'FontSize',16);
    set(ylab,'FontWeight','bold');

    set(hand,'FontSize',16);
    set(hand,'LineWidth',3);

    grid on;
    box on;
end



end










