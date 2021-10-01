function [closest,cval,err] = taft_findc(mat,val)

% [closest,cval,err] = taft_findc(mat,val)
%
% returns the index (closest), value (cval) and error (err) of the value closest to val
% in mat
% if val is an array or matrix, returns values for all elements of val
%
% Laurence Hunt
 
closest = nan(size(val));
cval    = nan(size(val));
err     = nan(size(val));
for i = 1:numel(val)
    [err(i),closest(i)] = min(abs(mat(:)-val(i)));
    cval(i)             = mat(closest(i));
end
