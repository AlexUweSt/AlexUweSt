a=TS{1,2};
b=TS{1,2};
k=0;


for lal=1:FileCount
a=TS{1,lal};
b=TS{1,lal};
k=0;
disp(lal);  
A = double.empty;
A1= double.empty;
problem=double.empty;
m=double.empty;
ss=50;
for i=1:size(a)-1
   b(i+1)=a(i)-a(i+1);
   if b(i+1)>=0 
       k=1;
   end
end


%vq = interp1(TS{1,1},R{1,1},TS{1,1});
plot(b)
j=1;
c=floor(size(a)/ss);
c=c(1);
for i=0:c-2
    suma=0;
    for q=1:ss
        suma=suma+b(i*ss+q+10);
    end
        
        
        
   if suma<(-0.009*ss) %0.009
    %disp(i*5);
    problem{j}=i*ss;
    j=j+1;
   end
end
problem1=cell2mat(problem);




for i=0:c-1
    if any(problem1(:) == i*ss);
        appnv=a(i*ss+1):(0.004*(1+rand)):a(i*ss+ss); %step 0.01
        A=vertcat(A,appnv');
    else 
        A=vertcat(A,a(ss*i+1:ss*i+ss));
        
    end
    
    
end


a=A;
for i=1:size(a)-1
   m(i+1)=a(i)-a(i+1);
end
enum=21;
c=floor(size(A)/enum);
c=c(1);
for i=0:c-1
   if any(a(i+1:i+enum)<-0.02)
       
      appnv=a(i*enum+1):(0.0035*(1+rand)):a(i*enum+enum);
        A1=vertcat(A,appnv');
    else 
        A1=vertcat(A,a(enum*i+1:enum*i+enum));    
    
    
    
end

end





%figure(2);
%a=A1;
%for i=1:size(a)-1
  % m(i+1)=a(i)-a(i+1);
%end
%plot(m)

 K=TS{1,lal};
R{1,lal} = interp1(K,R{1,lal},A1);  
  TS{1,lal}=A1;  
 
end



