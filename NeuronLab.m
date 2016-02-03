% adapted from Ryan Siciliano, 2012


%%% BUILD ONE PIECE AT A TIME: membrane, Na, K, n, KM
 am = @(v) 0.1*(v+35)./(1-exp(-(v+35)/10));
 bm = @(v) 4*exp(-0.0556*(v+60));
 an = @(v) 0.01*(v+50)./(1-exp(-(v+50)/10));
 bn = @(v) 0.125*exp(-(v+60)/80);
 m_inf = @(v) am(v)./(am(v) + bm(v));
 tau_m = @(v) 1./(am(v) + bm(v));
 n_inf = @(v) an(v)./(an(v) + bn(v));
 tau_n = @(v) 1./(an(v) + bn(v));
 ah = @(v) 0.07*(exp(-0.05*(v+60)));
 bh = @(v) 1./(1+exp(-0.1*(v+30)));
 h_inf = @(v) ah(v)./(ah(v) + bh(v));
 tau_h = @(v) 1./(ah(v) + bh(v));
 akm = @(v) 0.02./(1+exp(0.2*(-v-20)));
 bkm = @(v) 0.01*exp((-v-43)/18);
 km_inf = @(v) akm(v)./(akm(v) + bkm(v));
 tau_km = @(v) 1./(akm(v) + bkm(v));
 
Cm=0.01; % Membrane Capcitance uF/cm^2
dt=0.01; % ms Time Step ms
tmax = 500; % Duration of experiment
t=0:dt:tmax; % ms Time Array

El=-49.42; % mV Leakage reversal potential
ENa=55.17;
gbarNa=1.2;
gbarl=0.003;% mS/cm^2 Leakage conductance
EK=-72.14;
gbarK=0.36;
gbarkm=0.05;
V = zeros(size(t));
% I = zeros(size(t)); % External current applied per unit area (mA/cm^2)
m = zeros(size(t));
n = zeros(size(t));
h = zeros(size(t));
km= zeros(size(t));
% step_indices = find(t>50 & t<=100);
%  I(step_indices) = 0.2;
%  P = 20;
I = zeros(size(t)) + 0.14;
% ramp_start = .2; 
% ramp_end = 0.08;
% ramp_length = ramp_end - ramp_start; 
% ramp_steps = length(t)/2;
% I = zeros(size(t))+ramp_end;
% I(1:end/2) = ramp_start:ramp_length/(ramp_steps-1):ramp_end;
% pulse_indices = find(mod(t, P)<2);
% I(pulse_indices) = 0.07;
V(1)=-60; % Initial Membrane voltage
am(-80 : 0.01 : 80);
bm(-80 : 0.01 : 80);
m_inf(-80 : 0.01 : 80);
tau_m(-80 : 0.01 : 80);
m(1) = m_inf(V(1));
an(-80 : 0.01 : 80);
bn(-80 : 0.01 : 80);
n_inf(-80 : 0.01 : 80);
tau_n(-80 : 0.01 : 80);
n(1) = n_inf(V(1));
h_inf(-80 : 0.01 : 80);
tau_h(-80 : 0.01 : 80);
h(1) = h_inf(V(1));
akm(-80 : 0.01 : 80);
bkm(-80 : 0.01 : 80);
km_inf(-80 : 0.01 : 80);
tau_km(-80 : 0.01 : 80);
km(1) = km_inf(V(1));
for i=1:length(t)-1

 gl=gbarl;
 Il=gl*(V(i)-El);
gNa=gbarNa*(m(i)^3)*h(i);
INa=gNa*(V(i)-ENa);
gK=gbarK*n(i)^4;
IK=gK*(V(i)-EK);
gKM=gbarkm*km(i);
IKM=gKM*(V(i)-EK);

 %Euler method to find the next voltage value
 dVdt = (1/Cm) * (I(i) - (Il+INa+IK+IKM));
 V(i+1) = V(i) + dt*dVdt;
dmdt = (m_inf(V(i))-m(i))/tau_m(V(i));
m(i+1)= m(i)+dt*dmdt;
dndt = (n_inf(V(i))-n(i))/tau_n(V(i));
n(i+1)= n(i)+dt*dndt;
dhdt = (h_inf(V(i))-h(i))/tau_h(V(i));
h(i+1)= h(i)+dt*dhdt;
dkmdt = (km_inf(V(i))-km(i))/tau_km(V(i));
km(i+1)= km(i)+dt*dkmdt;

end
V_list = -80 : 0.01 : 80;
figure()
subplot(3,1,1)
plot(t, I)
ylabel('I(t)')
ylim([-0.2, 0.2])
subplot(3,1,2)
plot(t, V)
ylabel('V')
ylim([-100, 100])
xlabel('t (ms)')
subplot(3,1,3)
plot(t,m,'r')
hold on
plot(t,n,'g')
hold on 
plot(t,h,'y')
% figure();
% plot(V_list, am(V_list), 'r');
% hold on;
% plot(V_list, bm(V_list), 'b');
% figure();
% 
% plot(V_list, m_inf(V_list),'r')
% hold on
% plot(V_list,tau_m(V_list),'b');
% 
%  