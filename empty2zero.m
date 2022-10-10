function out = empty2zero(X)
    out = min([ max([X,0]),X ]); % - Convert "X = []" to "X = 0" while Keeping "X >< 0";
end