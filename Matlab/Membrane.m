clear all

fs = 44100;
k = 1 / fs;

lengthSound = fs;
drawThings = false;
drawSpeed = 100;

%% Membrane parameters
Lx = 0.3;     % Width
Ly = 0.3;     % Height
rho = 10; % Density
H = 0.001;  % Thickness
T = 50;      % Tension
E = 2e3;   % Young's modulus
nu = 0.3;   % Poissons ratio
D = E * H^3 / (12 * (1 - nu^2)); % will be divided by rho H in "createMembrane"

s0 = 0.1;
s1 = 0.01;
% [B, C, N, Nx, Ny, h, kappa, Dxx] = newCreatePlate (Lx, Ly, rho, H, D, s0, s1, k);
[N, Nx, Ny, h, c, kappa] = createMembrane (Lx, Ly, rho, H, T, D, s0, s1, k);
Nx = Nx - 1;
Ny = Ny - 1;

uNext = zeros (Nx, Ny);
u = zeros (Nx, Ny);
% u (0.5 * Nx + 0.5 * (Nx * Ny)) = 1;


exciterPosX = 0.3;
exciterPosY = 0.5;
%     rcW = floor(min(Nx,Ny) / 5);
rcW = floor(Nx / 3);
excitationMat = zeros(rcW+1, rcW+1);
scaler = (1-cos(2*pi*(0:rcW)/rcW)) * 0.5;
for x = 1:rcW+1
    excitationMat(x,:) = scaler(x) * (1-cos(2*pi*(0:rcW)/rcW)) * 0.5;
end

startIdxX = floor(Nx * exciterPosX - rcW/2) - 1;
startIdxY = floor(Ny * exciterPosY - rcW/2) - 1;
u(startIdxX + (0:rcW), startIdxY + (0:rcW)) = excitationMat;
% for i = 1 : rcW 
%    u((startIdxY + i) * Nx + startIdxX :(startIdxY + i) * Nx + startIdxX + rcW) = ...  
%         excitationMat(i,:);
% end  

uPrev = u;

vec = 3:Nx-2;
eVec = 2:Nx-1;
kinEnergy = zeros(lengthSound, 1);
potEnergyPlate = zeros(lengthSound, 1);
potEnergy1 = zeros(lengthSound, 1);
potEnergy2 = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);
out = zeros(lengthSound, 1);
lambdaSq = c^2 * k^2 / h^2;
muSq = kappa^2 * k^2 / h^4;

startT = 0;
Tinit = T;
for n = 1 : lengthSound
%     if (n > startT)
%         T = Tinit * (1 - 0.5 * ((n - startT)/lengthSound));
%     else
%         T = Tinit;
%     end
%     
    lambdaSq = T * k^2 / (rho * H * h^2);
    uNext(vec, vec) = (2 * u(vec,vec) - uPrev(vec, vec) + lambdaSq * (u(vec+1, vec) + u(vec-1, vec) + u(vec, vec+1) + u(vec, vec-1) - 4 * u(vec, vec))...
        - muSq * (u(vec+2, vec) + u(vec-2, vec) + u(vec, vec+2) +  u(vec, vec-2))...
        - 2 * muSq * (u(vec+1, vec+1) + u(vec+1, vec-1) + u(vec-1, vec+1) +  u(vec-1, vec-1))...
        + 8 * muSq * (u(vec+1, vec) + u(vec-1, vec) + u(vec, vec+1) +  u(vec, vec-1))...
        - 20 * muSq * u(vec,vec) + s0 * k * uPrev(vec,vec)...
        + 2 * s1 * k / h^2 * (u(vec+1, vec) + u(vec-1, vec) + u(vec, vec+1) + u(vec, vec-1) - 4 * u(vec, vec)...
        -(uPrev(vec+1, vec) + uPrev(vec-1, vec) + uPrev(vec, vec+1) + uPrev(vec, vec-1) - 4 * uPrev(vec, vec)))) / (1 + s0*k);
    
%     kinEnergy(n) = sum(sum(rho * H / 2 * h^2 * (1 / k * (u(eVec, eVec) - uPrev(eVec, eVec))).^2));
%     potEnergy1(n) = T / 2 * h^2 * sum(sum(1/h^2 * (u(eVec+1, eVec) - u(eVec, eVec)) .* (uPrev(eVec+1, eVec) - uPrev(eVec, eVec))...
%         + 1/h^2 * (u(eVec, eVec+1) - u(eVec, eVec)) .* (uPrev(eVec, eVec+1) - uPrev(eVec, eVec))));
%     potEnergy2(n) = D / (2 * h^2) * sum(sum((u(eVec+1, eVec) + u(eVec-1, eVec) + u(eVec, eVec+1) + u(eVec, eVec-1) - 4 * u(eVec, eVec))...
%         .* (uPrev(eVec+1, eVec) + uPrev(eVec-1, eVec) + uPrev(eVec, eVec+1) + uPrev(eVec, eVec-1) - 4 * uPrev(eVec, eVec))));
%     potEnergy2(n) = sum(sum(D / 2 * h^2 * 1/h^4 * (u(vec+1, vec) + u(vec-1, vec) + u(vec, vec+1) + u(vec, vec-1) - 4 * u(vec, vec))...
%         .* (uPrev(vec+1, vec) + uPrev(vec-1, vec) + uPrev(vec, vec+1) + uPrev(vec, vec-1) - 4 * uPrev(vec, vec))));
%     uNext = B * u + C * uPrev; 
%     uReshaped = reshape(reshape(u, Nx, Ny)', Nx * Ny, 1);
%     uPrevReshaped = reshape(reshape(uPrev, Nx, Ny)', Nx * Ny, 1);
%     kinEnergy(n) = sum(rho * H / 2 * h^2 * (1 / k * (u - uPrev)).^2);
% %     kinEnergy(n) = ((rho * H) / 2) * h^2 * sum(sum(1/k^2 * (u - uPrev).^2));
% %     potEnergy1(n) = sum(T / 2 * h^2 * 1/h^2 * (Dxy * u).^2 - 2 * Dxy * u);
% %     potEnergy2(n) = sum(D / 2 * h^2 * 1/h^4 * Dxx * (u - uPrev));
% %     potEnergy2(n) = D / (2 * h^2) * sum((Dxx * u) .* (Dxx * uPrev));
%     kinEnergyPlate(n) = ((rho * H) / 2) * h^2 * sum(sum(1/k^2 * (u - uPrev).^2));
% %     potEnergyPlate1(n) = T / 2 * h^2 * 1/h^2 * sum(((Dxy * u).^2 - 2 * Dxy * u) .* ((Dxy * uPrev).^2 - 2 * Dxy * uPrev));
%     potEnergyPlate1(n) = T / 2 * h^2 * 1/h^2 * sum(((u(2:end) - u(1:end-1)) .* (uPrev(2:end) - uPrev(1:end-1)))...
%         + ((uReshaped(2:end) - uReshaped(1:end-1)) .* (uPrevReshaped(2:end) - uPrevReshaped(1:end-1))));
%     potEnergyPlate(n) = D / (2 * h^2) * sum((Dxx * u) .* (Dxx * uPrev));
% 
    totEnergy(n) = kinEnergy(n) + potEnergy1(n) + potEnergy2(n);% + potEnergyPlate(n);% + potEnergy2(n);
    
    out(n) = u(floor(0.7 * Nx), floor(0.4 * Ny));
    if drawThings && mod(n,drawSpeed) == 0
%         subplot(3,1,1);
%         cla
%         hold on;
%         plot(kinEnergy(10:n))
%         plot(potEnergy1(10:n))
%         plot(potEnergy2(10:n))
%         subplot(3,1,2)
%         plot(totEnergy(10:n)/totEnergy(10)-1);% / totEnergy(10) - 1)
%         
%     %     zerosMat(2:Nx+1, 2:Ny+1) = reshape(uReshaped, Nx, Ny);
%     %      subplot(4,1,[3, 4])
%         subplot(3,1,3)
        mesh (u);
        zlim([-max(max(excitationMat)) * 0.5, max(max(excitationMat)) * 0.5])
        drawnow;
    end
    uPrev = u;
    u = uNext;
end
% if ~drawThings
    plot(out)
% end