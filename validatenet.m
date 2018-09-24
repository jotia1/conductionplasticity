function [  ] = validatenet( net )
%VALIDATENET Check a network definition is valid
%   Check a network definition is valid or raise an error if not.

assert(sum(net.group_sizes) == net.N, 'Group sizes do not sum to N');

% TODO: validate delays, variance and w have same connectivity


end

