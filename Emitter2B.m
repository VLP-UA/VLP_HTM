function [ B ] = Emitter2B( Emitters )
%Emitter2B    Compute B matrix from emitters data
%   Detailed explanation goes here

% B: matrix with emitter data
% $b_{ij} = \frac{m_j+1}{2\pi} ~ \delta_{ij}$ 

% Bv is a vector with the b_ij elements
Bv = ([Emitters(:).m] + 1)/(2*pi);

nv = numel(Bv);
B = sparse(1:nv,1:nv,Bv);
% B = diag(Bv);


end

