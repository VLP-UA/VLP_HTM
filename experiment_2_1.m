%% Computation of illumination by a set of light sources





%% Compute A and B matrices

% A: matrix with receiver data

% a_{ij} = A_{r,i} T_{s,i} \frac{n_i^2}{\sin^2(\Psi_{c,i})}
Av = [Receivers(:).Ar].*[Receivers(:).Ts].*([Receivers(:).n].*sin([Receivers(:).Psi]).^2);

% Compute the diagonal matrix with receiver data. 
%
% When there are many receivers, the direct use of diag(Av) generates an
% excessively large matrix:
% A = diag(Av);
% Better use a sparse matrix:
nv = numel(Av);
% Create the diagonal matrix using a sparse matrix:
A = sparse(1:nv,1:nv,Av);


% B: matrix with emitter data

% $b_{ij} = \frac{m_j+1}{2\pi} ~ \delta_{ij}$ 
Bv = ([Emitters(:).m] + 1)/(2*pi);

B = diag(Bv);

%% Compute F matrix

% Extract the position and orientation vectors from the HTMs

experiment_2_1_aux

Pr = A*F*B*[Emitters.Pt]';
vals = reshape(Pr,[nx ny]);

