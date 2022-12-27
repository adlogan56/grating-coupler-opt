function [ obj, obj_parts] = gratingObjectiveGeneral( gr_spec, S, L, center, sigma )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

z = zeros(1,length(S));
for iz = 2:length(S)
    z(iz) = z(iz - 1) + L(iz - 1) + gr_spec{1}.d(S(iz - 1));
end

obj_parts = zeros(1,length(gr_spec));

for ii = 1:length(obj_parts)
    %normalized gaussian
    gaussian = exp(-(z - center(ii)).^2/(2*sigma(ii)^2));
    gaussian = (1/sum(abs(gaussian).^2))*gaussian;
    
    [totalScatter, totalRef, E_field] = serialScattererFunctionGrating3(gr_spec{ii}.r,gr_spec{ii}.s,gr_spec{ii}.t,gr_spec{ii}.neff,gr_spec{ii}.lambda,S,L);
    
    obj_parts(ii) = abs(sum(totalScatter.*gaussian))^2;
end

obj = prod(obj_parts);

if length(obj_parts)>1
    obj = (obj*min(obj_parts))^(1/(length(obj_parts)+1));
end
    
    
end