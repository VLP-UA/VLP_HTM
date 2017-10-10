function [ Emitters ] = CreateEmittersArray2( n_Emitters, Pt, m, Lx, Ly, Z )
%CREATEEMITTERSARRAY2 Creates a regularly spaced 2 dimensional array of n_Emitters
%   Emitters = CreateEmittersArray( n_Emitters, Pt, m, Lx, Ly, Z )
%
%   n_Emitters  - number of emitters
%   Pt          - Transmited power in each emitter
%   m           - Lambertian order of the emitter
%   Lx, Ly      - Total length of the array in x and y directions, resp.
%   Z           - z coordinate (if real) or base HTM (HTM of the 1st emitter)
%
%   The emitters are place in a nx x ny array, where nx*ny= n_Emitters at
%   height Z, looking down (i.e., with z axis pointing downwards) when Z is
%   a real value or placed starting at the location defined by the HTM Z.

[ nx , ny ] = factorise2(n_Emitters);

Emitters = CreateEmittersArray(nx, ny, Pt, m, Lx, Ly, Z);