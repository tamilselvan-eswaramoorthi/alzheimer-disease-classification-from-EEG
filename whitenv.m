function [newVectors, whiteningMatrix, dewhiteningMatrix] = whitenv ...
    (vectors, E, D, s_verbose);

% Default value for 'verbose'
if nargin < 4, s_verbose = 'on'; end

% Check the optional parameter verbose;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end
if any (diag (D) < 0),
  error (sprintf (['[ %d ] negative eigenvalues computed from the' ...
		   ' covariance matrix.\nThese are due to rounding' ...
		   ' errors in Matlab (the correct eigenvalues are\n' ...
		   'probably very small).\nTo correct the situation,' ...
		   ' please reduce the number of dimensions in the' ...
		   ' data\nby using the ''lastEig'' argument in' ...
		   ' function FASTICA, or ''Reduce dim.'' button\nin' ...
		   ' the graphical user interface.'], ...
		  sum (diag (D) < 0)));
end

whiteningMatrix = inv (sqrt (D)) * E';
dewhiteningMatrix = E * sqrt (D);

if b_verbose, 
%     fprintf ('Whitening...\n');
end
newVectors =  whiteningMatrix * vectors;

if ~isreal(newVectors)
  error ('Whitened vectors have imaginary values.');
end

if b_verbose
%   fprintf ('Check: covariance differs from identity by [ %g ].\n', ...
%     max (max (abs (cov (newVectors', 1) - eye (size (newVectors, 1))))));
end
