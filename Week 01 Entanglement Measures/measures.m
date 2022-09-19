% defining the states
psi = [1; 0; 0; 1] / sqrt(2);
rho = zeros(4,4);
rho(1,1) = .5;
rho(4,4) = .5;

psi = psi * transpose(psi);

fprintf("EE(psi) = %.5f\n", EntanglementEntropy(psi));
fprintf("EE(rho) = %.5f\n", EntanglementEntropy(rho));
fprintf("EoF(psi) = %.5f\n", EoF(psi));
fprintf("EoF(rho) = %.5f\n", EoF(rho));
fprintf("LN(psi) = %.5f\n", LogratihmicNegativity(psi));
fprintf("LN(rho) = %.5f\n", LogratihmicNegativity(rho));

% entanglement entropy
function S = EntanglementEntropy(rho)
  TrB = PartialTraceB(rho);
  ev = eigs(TrB);

  S = -1 * (ev(1) * log2(ev(1)) + ev(2) * log2(ev(2)));
end

% Entanglement of formation
function E = EoF(rho)
  E = shannon (1/2 + 1/2 * sqrt(1 - concurrence(rho)^2));
end

function LN = LogratihmicNegativity(rho)
  [mat, res] = sqrtm(ctranspose(PartialTransposeA(rho)) * PartialTransposeA(rho));
  LN = log2(trace(mat ));
end

% helper functions
% stuff like classical entropy, partial transpose, partial trace

% Shannon entropy
function h = shannon(x)

  if (x == 0 || x == 1)
    h = 0;
  else
    h = -1 * (x .* log2(x) + (1-x) .* log2(1-x));
  end
end

% concurrence
function C = concurrence(rho)
  sigm_y = [0, -1i;1i, 0];
  sigm_y_2 = kron(sigm_y, sigm_y);
  rho_tilde = sigm_y_2 * conj(rho) * sigm_y_2;

  [sqrt_rho, residual] = sqrtm(rho);
  [R, residual] =  sqrtm(sqrt_rho * rho_tilde * sqrt_rho);

  l = eigs(R);
  C = max(0, l(1) - l(2) - l(3) - l(4));
end

% function taking the trace of B for a Operator in the Space AxB
function rho_A = PartialTraceB(rho)
  rho_A = zeros(2,2);

  rho_A(1,1) = rho(1,1) + rho(2,2);
  rho_A(2,2) = rho(3,3) + rho(4,4);

  rho_A(1,2) = rho(1,3) + rho(2,4);
  rho_A(2,1) = rho(3,1) + rho(4,2);
end

% function taking the partial transpose over the first system
function rho_transp = PartialTransposeA(rho)
  rho_transp = rho;
  rho_transp(1:2, 3:4) = rho(3:4, 1:2);
  rho_transp(3:4, 1:2) = rho(1:2, 3:4);
end

