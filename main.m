clear;
tic
% Load mesh and plot:
load('Mesh.mat');
[TR, Db, nC, d, nE, dNodes, fNodes, s, m, vol_T, mp_T, r_c4n, theta_c4n] = build_mesh(c4n, n4e, R);

N = 20; % size of K_n. Must be even

% Indices alpha and beta of K range from -N/2 to N/2
BETA  = -N/2:N/2;
ALPHA = -N/2:N/2;

% lattice in complex plane that contains spectral parameter:
M = 20;
min_real = 1;
max_real = 3;
min_imag = -0.1;
max_imag = 0.01;
adjustment = (max_real-min_real)/(max_imag-min_imag);
X = linspace(min_real, max_real, adjustment*M);
Y = linspace(min_imag, max_imag, M);
[XX,YY] = meshgrid(X,Y);
L = XX + 1i*YY;

% Square Root of N-Matrix and projection P_0:
p0 = 1; % 0th entry of P_0
NN_diag = zeros(N+1,1);
for l=(-N/2:N/2)
   NN_diag(l+N/2+1) = sqrt(abs(l)) ;
end
NN_diag(N/2+1) = p0;
NN_inv         = diag(1./NN_diag);
NN = diag(abs(-N/2:N/2));
P0 = zeros(N+1);
P0(N/2+1,N/2+1) = p0;

% Initialize Determinant
Determinant = zeros(size(L));

% FEM Dirichlet boundary data:
tu_D = u_D(r_c4n,theta_c4n,ALPHA,r_outer,R);
phi  = u_D(r_c4n,theta_c4n,BETA,r_outer,R);

%% Main Loop:

disp('Computing...')
parfor i=1:size(L,1)*size(L,2) % loop over all z in L
    z = L(i); 

    % Definition of H:
    H_diag = -z*besselh(abs(ALPHA)-1,z*R)./besselh(abs(ALPHA),z*R);
    HH     = diag(H_diag);
        
    % Finite element approximation for M_inner:
    u           = zeros(nC,N+1);
    S           = s-z^2*m;   % weak version of -Δ-z^2
    b           = -S*tu_D;
    u(fNodes,:) = S(fNodes,fNodes)\b(fNodes,:); 
    u           = u+tu_D;

    M_inner        = phi'*S*u;
    A = 0.5*(eye(N+1) - P0 + R*NN_inv*(HH + M_inner)*NN_inv);
    Determinant(i) = det(A);
end
disp('Done!')
toc
%% Plot Results in complex plane:

figure
% Contour plot of log(|det(I+H+J+K)|):
[Mat c] = contour(real(L), imag(L), log(abs(Determinant)), 70);
colorbar;
legend('log(|det(I+H+J+K)|)','Location','southwest');

% Resonances are approximated by local minima of Determinant:
Resonances = islocalmin(abs(Determinant),1) & islocalmin(abs(Determinant),2);
resonance_values = L(Resonances)

% Dirichlet EVs of the closed chamber (only relevant for circular resonator):
EVs = [2.8531, 1.3360, 2.1287];
EVs = cast(EVs,'like',L);
hold on 
plot(EVs,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'MarkerSize',7, 'DisplayName', 'D-EVs')
plot(L(Resonances),'o','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0 0 0],'MarkerSize',8, 'DisplayName', 'Resonances')
hold off

