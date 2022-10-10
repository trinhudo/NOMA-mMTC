function out = empty2Inf(X)
    out = max([ min([X,Inf]),X ]); % - Convert "X = []" to "X = Inf" while Keeping "X >< 0";
end