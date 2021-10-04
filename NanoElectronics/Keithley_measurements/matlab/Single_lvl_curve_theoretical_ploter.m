%SLM.current(1,1,1,1);
figure(1)
x=-0.8:0.0001:0.8;%0.0001

h=zeros(1601);
t=0.24
g=0.005
t=0.245
g=0.005
h=semilogy(x,abs(SLM.current(t,g/2,g/2,x)),'LineWidth',2,'Color','black');
alpha(0.5)
pbaspect([1 1 1])
xlim([-0.8 0.8])
ylim([1e-10 1e-6])
%%
%hold on
%h=zeros(1601);
%t=0.24
%g=0.005
%t=0.24
%g=0.005
%h=semilogy(x,abs(SLM.current(t,g/2,g/2,x)),'LineWidth',2,'Color','red');
%alpha(0.5)
%pbaspect([1 1 1])
%xlim([-0.8 0.8])
%ylim([1e-10 1e-6])

%%
%saveas(gcf,'C:/Users/Filip/Desktop/myGray.png')

figure(2)
x=-0.8:0.0001:0.8;

h=zeros(1601);
t=0.68
g=0.039
t=(0.6+0.81)/2
g=(0.028+0.046)/2
h=semilogy(x,abs(SLM.current(t,g/2,g/2,x)),'LineWidth',2,'Color','black');
alpha(0.5)
pbaspect([1 1 1])
xlim([-0.8 0.8])
ylim([1e-10 1e-6])
%%

figure(3)
x=-0.8:0.0001:0.8;

h=zeros(1601);
t=0.96
g=0.025
t=(0.78+1.18)/2
g=(0.014+0.030)/2
h=semilogy(x,abs(SLM.current(t,g/2,g/2,x)),'LineWidth',2,'Color','black');
alpha(0.5)
pbaspect([1 1 1])
xlim([-0.8 0.8])
ylim([1e-10 1e-6])

