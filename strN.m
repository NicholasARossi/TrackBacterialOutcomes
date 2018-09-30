% This function takes a number x and returns a string of length N with the
% appropriate number of zeros in front of it.
% e.g. strN(7, 3) returns '007'

function c = strN(x, N)

L = length(num2str(x));
if L > N
    disp('Warning: strN has length of string larger than N');
end
c = num2str([zeros(1,N-L) x]);
c(c==' ') = '';                     % cut out whitespace
