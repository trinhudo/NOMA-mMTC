% [1;    [1*ones(1,N);
%  2; --> 2*ones(1,N);
%  3]     3*ones(1,N)];
function out = vec2mat(vec,N)
    out = cell2mat( arrayfun( @(x) vec(x)*ones(1,N),(1:length(vec))','UniformOutput',false ) );
end