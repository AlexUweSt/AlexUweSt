


for i=1:FileCount
  n=size(TS{1,i});
p=TS{1,i}(n(1));
v=p(1);
inc=v/n(1);
A=TS{1,i}(1):inc/8:TS{1,i}(n);
[TS{1,i}, index] = unique(TS{1,i});
R{1,i} = interp1(TS{1,i},R{1,i}(index),A)';  
  TS{1,i}=A';  
end