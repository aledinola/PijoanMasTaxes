function [mystruct] = pack_into_struct(varargin)
% pack_into_struct packs input arguments into a structure
% EXAMPLE
% a = 1; b= 2; c = [3,4]; d= 'pippo';
% par = pack_into_struct(a,b,c);
% will result in a structure par with fields a,b,c and d.

n = numel(varargin);

for ii = 1:n
    var_name  = inputname(ii);
    var_num = varargin{ii};
    mystruct.(var_name) = var_num;
end

end %end function

