%Experiment 2_2

% Computation of illumination by a set of light sources

% Same operations as experiment_2_1, but using Receiver2A, Emitter2B and
% EmitRec2F functions.

A = Receiver2A(Receivers);

B = Emitter2B(Emitters);

F = EmitRec2F(Emitters, Receivers);


Pr = A*F*B*[Emitters.Pt]';

vals = reshape(Pr,[nx ny]);

