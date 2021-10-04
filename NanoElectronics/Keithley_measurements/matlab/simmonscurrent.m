 %% parameters
 
hold on
beta=1; %junction shape parameter
m = 9.10938291E-31; %electron mass [kg]
hbar = 1.054571726E-34; %Planck constant [Js]
e=1.60217662E-19; % electron charge in Coulombs [C]
fi=5.1;
fie=fi*e; %potential [J]
V=0.1; %applied voltage to curcuit [V]
r=1e-10; %[m]
Amp=6.241E18; %Coulomb
G0=7.7480E-5; %quantum of conductance [S]
R0=1/G0; %quantum of resistance [Ohm]
R=2000; %pre-resistor[Ohm]
Rall=R+R0; %resistance of resistor + one atomic contact 
K=2.38; %relative dielectric index of toluene

%% basic simmons model without image forces and dielectric index, and low voltage regime
%J=V*(gama*sqrt(fi)/d)*exp(-A*d*sqrt(fi))

A=2*beta*sqrt(2*m/hbar^2);
gama=(e*sqrt(2*m)/(4*beta*pi^2*hbar^2));
coef=@(d)gama*(sqrt(fie)/d)*exp(-A*d*sqrt(fie))/Amp; % %J=V*coef(d), plus added division by Amp to convert to amperes

%j2=(V-Rall*j1)*coef , equation for current; 
% V-Rall*j1 represents voltage drop between the tips

j1=@(t)V/(1/coef(t)+Rall); %current through the junction
Cond=@(q)1/(((V/j1(q)-2000))/R0); % conversion to conductance of the junction including the monoatomic contact
%% plot of the basic simmons model

x=0:1e-11:2.5e-9; % array of point for x-axis
%y=arrayfun(j1,x); 
i=arrayfun(Cond,x); %applying the array to the 'Cond' function

for m= 0:1e-11:2.5e-9 %detection of zero position on histogram
    if Cond(m)<0.99 
        x0=m;
        break;
    end    
end

%x0=0; %obrisi ovo
%figure(11);

%plot(x-x0,log10(i)); %for normal plotting without semilog scale 
%semilogy(x-x0,i,'blue','LineWidth',1);%x-1.65e-9 %for ploting on semilog
%scale
%% simmons model with image forces and dielectric constant K


s1=6/(K*fi); %formula No.49 in article(original)
s2=@(d)(d*1e10-s1);
deltas=@(d)(s2(d)-s1);
%modified potential; distance 'd' in angstrems,fi potential in eV
fi2=@(d)(fi-((5.75)/(K*(s1-s2(d))))*log(s2(d)*(d*1e10-s1)/(s1*(d*1e10-s2(d)))));
%fi2 potential in eV
coef2=@(d)3.16e10*(sqrt(fi2(d))*(1/deltas(d)))*exp(-1.025*deltas(d)*sqrt(fi2(d))); 
%(1e10/deltas(d)) converts from units of angstrems to [m]
j2=@(t)V/(1/coef2(t)+Rall);
Cond2=@(q)1/(((V/j2(q)-2000))/R0);

%% plot of the modified simmons model


x2=0:1e-11:2.5e-9;
i2=arrayfun(Cond2,x2);

for m= 0:1e-11:2.5e-9 %detection of zero position on histogram
    if Cond2(m)<0.9 
        x0=m;
        break;
    end    
end

figure(1);
%plot(x2-x0,log10(i2),'red');
semilogy(x2-x0,i2,'red','LineWidth',1); %for use with semilog scale

