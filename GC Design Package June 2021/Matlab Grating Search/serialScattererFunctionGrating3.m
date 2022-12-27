function [totalScatter, totalRef, E_field] = serialScattererFunctionGrating3(R,S,T,n_eff,lambda,scatterers,L)

L_eff  = [L 0];     % Add a placeholder 0 length after last real scatterer

N = length(scatterers);
k = 2*pi*n_eff./lambda;

E_field = zeros(2,2*N+1);    % First row is forward-propagating, second is back
E_field(1,2*N+1) = 1;

for ii = N:-1:1
    
    s_type = scatterers(ii);
    r = R(s_type);
    t = T(s_type);
    l = L_eff(ii);
    
    % Each step starts with the field just before a scatterer
    
    mat_prop = [exp(-1i*k*l) 0; 0 exp(1i*k*l)];      % Extrapolate light after previous scatterer from propagation in wg
    mat_scatter = (1/t)*[1 -r; r (t^2 - r^2)];      % Extrapolate light before scatterer from light after scatterer
    
    E_field(:,(2*ii)) = mat_prop*E_field(:,(2*ii)+1);
    E_field(:,(2*ii)-1) = mat_scatter*E_field(:,(2*ii));
end

totalRef = E_field(2,1)/E_field(1,1);

totalScatter = zeros(1,N);

for ij = 1:N
    
    s_type = scatterers(ij);
    s = S(s_type);
    
    totalScatter(ij) = (s/E_field(1,1))*(E_field(1,(2*ij)-1) + E_field(2,(2*ij)));
end

end



