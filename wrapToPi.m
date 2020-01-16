%% Copyright (C) 2015-2020 Oscar Monerris Belda
%% 
%% This program is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*- 
%% @deftypefn {Function File} {@var{xwrap} =} wrapToPi (@var{x})
%%
%% Wraps X into the [-pi to pi] interval
%%
%% @var{x}: angle(s) in radians (single value, vector or ND-matrix).  
%%
%% @var{xwrap}: output value(s) in the range [-pi .. pi] radians.
%% The interval [-pi .. pi] is a closed interval: values equal to
%% negative odd multiples of -pi are mapped to -pi, values equal to
%% an odd multiple of +pi are mapped to pi.
%%
%% Example:
%% @example
%% wrapToPi ([-3*pi, -pi, -pi-1, 0; pi-1, pi, pi+1, 3*pi])
%% ans =
%%  -3.14159  -3.14159   2.14159   0.00000
%%   2.14159   3.14159  -2.14159   3.14159
%% @end example
%%
%% @seealso{wrapTo180, wrapTo360, wrapto2Pi}
%% @end deftypefn

function xwrap = wrapToPi(x)

  xwrap = rem (x, 2*pi);
  idx = find (abs (xwrap) > pi);
  xwrap(idx) = xwrap(idx)- 2*pi * sign (xwrap(idx));

end

%!test
%! x = -9:0.1:9;
%! xw = wrapToPi (x);
%! assert (sin (x), sin (xw), 8 * eps)
%! assert (cos (x), cos (xw), 8 * eps)
%! assert (! any (xw < -pi))
%! assert (! any (xw > pi))

%% Test Matlab compatibility as regards closed interval (esp. left)
%!test
%! assert (wrapToPi ([-3*pi, -pi, -pi-1, 0; pi-1, pi, pi+1, 3*pi]), ...
%!                   [-pi, -pi, pi-1, 0.00000 ; ...
%!                    pi-1, pi, -pi+1, pi], 1e-13)
