function L = semispectrum_1DSD(order,nbrCells,upwindPar,doplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Construct the linear operator (matrix) L of the Spectral Difference method for the 1D linear
%%  advection equation.
%%
%%  Inputs:
%%          
%%  -order of the scheme from 1 to 6.
%%   (1st-order SD coresponds to the classical FD method!!!
%%
%%  - upwind parameter
%%  - uniform or perturbed grid
%%
%%
%%  Author: Matteo Parsani
%% 
%%          parsani.matteo@gmail.com
%%
%% Date: 03-23-2011 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Preliminaries
% %%%%%%%%%%%%%%%
% clear all;
% close all;
% clc;


% Flux and solution points coordinate
%************************************
switch order
    case 1
        % 1st order
        flxPnts = [ -1.0 ; 1.0 ];
        solPnts = [ 0 ];
    case 2
        % 2nd order
        flxPnts = [ -1.0 ; 0.0 ; 1.0 ];
        solPnts = [ -1.0 ;       1.0 ];
    case 3
        % 3rd order
        flxPnts = [ -1.0 ; -0.58 ; 0.58 ; 1.0 ];
        solPnts = [ -1.0 ;         0.58 ; 1.0 ];
    case 4
        % 4th order
        flxPnts = [ -1.0 ; -0.78 ; 0.0 ; 0.78 ; 1.0 ];
        solPnts = [ -1.0 ; -0.78 ;       0.78 ; 1.0 ];
    case 5
        % 5th order
        flxPnts = [ -1.0 ; -0.83 ; -0.36 ; 0.36 ; 0.83 ; 1.0 ];
        solPnts = [ -1.0 ; -0.83 ;         0.36 ; 0.83 ; 1.0 ];
    case 6
        % 6th order
        flxPnts = [ -1.0 ; -0.88 ; -0.53 ; 0.0 ; 0.53 ; 0.88 ; 1.0 ];
        solPnts = [ -1.0 ; -0.88 ; -0.53 ;       0.53 ; 0.88 ; 1.0 ];
    otherwise
        disp('Error: min order 1, max order 6');
end


% Number of solution and flux points
%***********************************
nbrflxPnts = length(flxPnts);
nbrSolPnts = length(solPnts);


% Matrices for the calculation of the solution at the flux points 
%****************************************************************
extrSolToFlx = ones(nbrflxPnts,nbrSolPnts);

for iFlx = 1:nbrflxPnts
    for iSol = 1:nbrSolPnts
        for iFactor = 1:nbrSolPnts
            if iFactor ~= iSol
                extrSolToFlx(iFlx,iSol) = extrSolToFlx(iFlx,iSol)*(flxPnts(iFlx)-solPnts(iFactor))/(solPnts(iSol)-solPnts(iFactor));
            end
        end
    end
end

%% NOTE: Multipliyng extrSolToFlx(iFlx,iSol) with sol(cellID,iSol) we
%% obtain the contribution of the iSol solution to the solution at the iFlx
%% point.



% Riemann flux parameters. 
% If fully upwind (upwindPar =1)
facLeft = 0.5*(1.0 + upwindPar); % Contribution of the left solution. 
facRigh = 0.5*(1.0 - upwindPar); % Contribution of the right solution

compFluxCurrCell = zeros(nbrflxPnts,nbrSolPnts);
compFluxLeftCell = zeros(nbrflxPnts,nbrSolPnts);
compFluxRighCell = zeros(nbrflxPnts,nbrSolPnts);


% Flux at the internal flux points
for iFlx = 2:nbrflxPnts-1
    for iSol = 1:nbrSolPnts
        compFluxCurrCell(iFlx,iSol) = extrSolToFlx(iFlx,iSol);  % We assume convective velocity = 1, i.e. the flux function is 
                                                                % f= a*u
    end
end

%% NOTE: Multiplying compFluxCurrCell(iFlx,iSol) with  sol(cellID,iSol) we
%% obtain the contribution of the iSol solution to the flux at the iFlx
%% point. 
%% compFluxCurrCell(iFlx,iSol) = extrSolToFlx(iFlx,iSol) only for internal
%% flux points because the flux at the face flux points must be corrected
%% with the solution computed with the Riemann solver.


% Flux at the boundaries flux points
for iSol = 1:nbrSolPnts
    % left face
    compFluxLeftCell(1         ,iSol) = facLeft*extrSolToFlx(nbrflxPnts,iSol); 
    compFluxCurrCell(1         ,iSol) = facRigh*extrSolToFlx(1         ,iSol);
    
    % right face
    compFluxCurrCell(nbrflxPnts,iSol) = facLeft*extrSolToFlx(nbrflxPnts,iSol);
    compFluxRighCell(nbrflxPnts,iSol) = facRigh*extrSolToFlx(1         ,iSol);
    
    
end



% Matrix for the calculation of the flux derivative at the solution points
%*************************************************************************
derivFlxInsolPnts = zeros(nbrSolPnts,nbrflxPnts);

for iSol = 1:nbrSolPnts
    for iFlx = 1:nbrflxPnts
        derivFlxInsolPnts(iSol,iFlx) = 0.0;
        for iTerm = 1:nbrflxPnts
            if iTerm ~= iFlx
                term = 1.0/(flxPnts(iFlx)-flxPnts(iTerm));
                for iFactor = 1:nbrflxPnts
                    if (iFactor ~= iTerm) && (iFactor ~= iFlx)
                        term = term*(solPnts(iSol)-flxPnts(iFactor))/(flxPnts(iFlx)-flxPnts(iFactor));
                    end
                end
                derivFlxInsolPnts(iSol,iFlx) = derivFlxInsolPnts(iSol,iFlx) + term;
            end
        end
        % Factor 2 comes from Jacobian determinant
        derivFlxInsolPnts(iSol,iFlx) = derivFlxInsolPnts(iSol,iFlx)*2.0;
    end
end



% Compute block matrices
DMm1  = derivFlxInsolPnts*compFluxLeftCell;
DM0   = derivFlxInsolPnts*compFluxCurrCell;
DMp1  = derivFlxInsolPnts*compFluxRighCell;

[p,q] = size(DM0); %size of each small matrix



% Construct the homogeneous part of matrix L (u^prime = L * u) using blktridiag.m function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = full(blktridiag(DM0,DMm1,DMp1,nbrCells));


% Add boundary condition for the first cell.
% This term is nothing than DMm1 and DMp1 placed at the right position.
% DMm1 affects the 1st cell.
% DMp1 affects the last cell.
% For a pure upwind scheme, the last cell is not affected by DMp1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = length(L(:,1));

% Apply BC to the first cell --> DMm1
 for i = 1:p
     for j = 1:p
         L(i,tmp-p+j) = DMm1(i,j);
     end
 end
 
% Apply BC to the last cell --> DMp1
 for i = 1:p
     for j = 1:p
         L(tmp-p+j,j) = DMp1(i,j);
     end
 end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% % Create mesh
% %%%%%%%%%%%%%
% xMin = 0;
% xMax = 1;
% 
% cellXMin(1)         = xMin;
% cellXMax(nbrCells)    = xMax;
% 
% deltaX = (cellXMax(nbrCells) - cellXMin(1))/nbrCells;
% 
% for iCell = 2:nbrCells
%     cellXMin(iCell) = xMin + deltaX*(iCell - 1 + perturbFactor*random('unif',-1,1));
%     %$cellXMin(iCell) = xMin + deltaX*(iCell - 1);
%     cellXMax(iCell-1) = cellXMin(iCell);
% end
% 
% for iCell = 1:nbrCells
%     cellXCen(iCell) = 0.5*(cellXMax(iCell) + cellXMin(iCell));
%     dX(iCell)   = cellXMax(iCell) - cellXMin(iCell);
% end
%
%  % Multiply by cell's size.
%  %%%%%%%%%%%%%%%%%%%%%%%%%%
%  for iCell = 1: nbrCells
%      for jSolPnts = 1: nbrSolPnts
%          deltaXFactor(jSolPnts+(iCell-1)*nbrSolPnts) = 1./dX(iCell);
%      end
%  end
% 
%  
% L = bsxfun(@times,L,deltaXFactor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
SpectrumSD = eig(-L);
 

if doplot
    % Plot eigenvalues for checking their position in the complex plane
    %******************************************************************
    figure;
    plot(real(SpectrumSD(:)),imag(SpectrumSD(:)),'o','markersize',8, 'MarkerEdgeColor','k',...
                                                                 'MarkerFaceColor','k');

    h = gca;
    title_handle = title(['SD spectrum for the linear advection equation [order=(p+1)= ',int2str(order),']']);
    set(title_handle,'FontSize',16);
    set(title_handle,'FontWeight','bold','interpreter', 'latex');

    xlabel('Re(eigenvalue)');
    xlab = get(h, 'xlabel');
    set(xlab,'FontSize',16);
    set(xlab,'FontWeight','bold','interpreter', 'latex');

    ylabel('Im(eigenvalue)');
    ylab = get(h, 'ylabel');
    set(ylab,'FontSize',16);
    set(ylab,'FontWeight','bold','interpreter', 'latex');

    set(h,'FontSize',16);
    set(h,'LineWidth',3);

    grid on;
    box on;
end


% Pass the correct matrix to convex_form.m (u^' = L u) 
L = -L;
                                                   





end
