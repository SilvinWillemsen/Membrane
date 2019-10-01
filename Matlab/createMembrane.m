function [N, Nx, Ny, h, c, kappa] = createMembrane (Lx, Ly, rho, H, T, D, s0, s1, k)
%     [B, C, N, Nx, Ny, h, c, kappa, D, Dxy]
    kappa = sqrt(D / (rho * H));
    c = sqrt(T / (rho * H));
    LWRatio = Lx/Ly;
%     h = 2*sqrt(k*(s1^2+sqrt(kappa^2+s1^2))); % need to adapt to include tension
    h = 2 * sqrt ((c^2 * k^2 + 4 * s1 * k + sqrt((c^2 * k^2 + 4 * s1 * k)^2 + 4 * kappa^2 * k^2))/2);
    Nx = floor(Lx/h);
    Ny = floor(Ly/h);
    
    h = max(Lx/Nx, Ly/Ny);
    N = (Nx-1)*(Ny-1);
    
%     phi = 2 * s1 * k / h^2;
%     % generate scheme matrices
%     Dxx = sparse(toeplitz([-2;1;zeros(Nx-3,1)]));
%     Dyy = sparse(toeplitz([-2;1;zeros(Ny-3,1)]));
%     
%     Dx = sparse(1:Nx-1, 1:Nx-1, -1) + sparse(1:Nx-2, 2:Nx-1, ones(1, Nx-2), Nx-1, Nx-1);
%     Dy = sparse(1:Nx-1, 1:Nx-1, -1) + sparse(1:Nx-2, 2:Nx-1, ones(1, Nx-2), Nx-1, Nx-1);
%     
%     Dxy = kron(speye(Nx-1), Dy)+kron(Dx, speye(Ny-1)); 
%     D = kron(speye(Nx-1), Dyy)+kron(Dxx, speye(Ny-1)); 
%     DD = D*D/h^4; 
%     B = sparse((2*speye(N) + c^2*k^2*D/h^2 - kappa^2*k^2*DD - phi * D)/(1 + s0*k));
%     C = -((1-s0*k - phi * D)/(1+s0*k)).*speye(N); 
    
end