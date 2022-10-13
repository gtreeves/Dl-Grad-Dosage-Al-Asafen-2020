function v = strfindiDU(cellstr,pattern)
%Just like "strfind" but gives true-false (rather than indices)
%
%


[m,n] = size(cellstr);

v = false(m,n);
for i = 1:m
	for j = 1:n
		v(i,j) = strcmpi(cellstr{i,j},pattern);
	end
end