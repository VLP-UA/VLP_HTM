% VLP_HTM   
%
%   Matlab library for applying Homogeneous Transformation Matrices to
%   Visible Light Positioning.
%
%   Generation of A, B and F matrices:
%   Receiver2A            - Computes the A matrix from an array of receivers struct data
%   Emitter2B             - Compute B matrix from emitters data
%   EmitRec2F             - Compute the F matrix from emitter and receiver data.
%
%   Channel modelling:
%   H0_ER                 - Compute OWC channel gain
%
%   Sensor displacement:
%   vlpMoveSensor         - Displaces a sensor or a receiver array in 3D space
%
%   Data generation:
%   newEmitters           - Creates new Emitters
%   newReceivers          - Creates new Receivers
%   CreateEmittersArray   - Creates a regularly spaced 2 dimensional array of nx x ny emitters
%   CreateEmittersArray2  - Creates a regularly spaced 2 dimensional array of n_Emitters
%   CreateReceiversArray  - Creates a regularly spaced array of nx x ny receivers in a area of Lx x Ly
%   CreateReceiversArray2 - Creates a regularly spaced array of n_Receivers receivers in a area of Lx x Ly
%
%   Demos:
%   demo1                 - Demonstrates the usage of HTM based VLP code
%   DefineSpace           - Define variables for the simulation space
%   sample1               - Sample code to demonstrate H0_ER function usage
%   experiment1           - Plot of power received in the plane from a source using H0_ER
%   experiment2           - Plot radiation graph from Lambertian source using H0_ER
%   Exp_20170411_a        - Compute the received power in a set of receivers from a set of emitters using A, B and F matrices.
%
%   Miscelaneous:
%   factorise2            - Computes two factors that, multiplied, result in n
