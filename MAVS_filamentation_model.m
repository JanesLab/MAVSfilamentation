function MAVS_filamentation_model
%Simulates the filamentation and downstream interferon-stimulated gene 
%(ISG) induction for three different MAVS alleles subject to different
%extents of cleavage by enteroviral 3C proteinase.  In the model, the Glu93Ala271
%mutant is uncleavable, the Glu93Gln271 allele can be cleaved one way
%(MAVSclv1), and the Gln93Gln271 allele can be cleaved one of two ways
%(MAVSclv1 and MAVSclv2).

clear all
close all hidden

%Set all initial conditions to zero
MAVS_IC = 0;
polyMAVS_IC = 0;
ISGs_IC = 0;
C3pro_IC = 0;
MAVSclv1_IC = 0;
MAVSclv2_IC = 0;
Network_IC = [MAVS_IC polyMAVS_IC ISGs_IC C3pro_IC MAVSclv1_IC MAVSclv2_IC];

%Define time interval and model parameters
tspan = 250;    %arbitrary time units, which propagate to rate parameters
k_act = 1;      %"average" activation rate
k_poly = 1000;  %"fast" polymerization rate, consistent with intramembrane association
k_deg = 0.1;    %"slow" degradation rate
k_txn = 1;      %"average" transcription rate
k_clv = 0.1;    %"slow" proteolytic cleavage rate
phi = 1;        %default negative-feedback strength
gamma = 1;      %ratio of negative-feedback strengths for MAVSclv2/MAVSclv1
filament = 800; %estimate of MAVS filament size, based on four CARD domains of 
                %MAVS causing an axial rise of ~5 angstroms
                %(https://www.ncbi.nlm.nih.gov/pubmed/25018021) and a 
                %maximum filament size of 100 nm
                %(https://www.ncbi.nlm.nih.gov/pubmed/23273991)
                %100 nm/ 0.5 nm * 4 = 800 MAVS monomers
Model_Params = [k_act k_poly k_deg k_txn k_clv phi gamma filament];

%Simulate time courses
[Glu93Ala271_time,Glu93Ala271_Network_soln] = ode15s(@(t,y) NetworkODE(t,y, ...
    Model_Params,'Glu93Ala271'),[0 tspan],Network_IC,[]);
[Glu93Gln271_time,Glu93Gln271_Network_soln] = ode15s(@(t,y) NetworkODE(t,y, ...
    Model_Params,'Glu93Gln271'),[0 tspan],Network_IC,[]);
[Gln93Gln271_time,Gln93Gln271_Network_soln] = ode15s(@(t,y) NetworkODE(t,y, ...
    Model_Params,'Gln93Gln271'),[0 tspan],Network_IC,[]);

%Plot solutions
figure(1)
subplot(6,1,1)
semilogy(Glu93Ala271_time,Glu93Ala271_Network_soln(:,1),'r');
hold on
semilogy(Glu93Gln271_time,Glu93Gln271_Network_soln(:,1),'g');
semilogy(Gln93Gln271_time,Gln93Gln271_Network_soln(:,1),'b');
ylabel('MAVS');
axis([0 tspan 1e-3 1e3])
subplot(6,1,2)
semilogy(Glu93Ala271_time,Glu93Ala271_Network_soln(:,2),'r');
hold on
semilogy(Glu93Gln271_time,Glu93Gln271_Network_soln(:,2),'g');
semilogy(Gln93Gln271_time,Gln93Gln271_Network_soln(:,2),'b');
ylabel('polyMAVS');
axis([0 tspan 1e-3 1e3])
subplot(6,1,3)
semilogy(Glu93Ala271_time,Glu93Ala271_Network_soln(:,3),'r');
hold on
semilogy(Glu93Gln271_time,Glu93Gln271_Network_soln(:,3),'g');
semilogy(Gln93Gln271_time,Gln93Gln271_Network_soln(:,3),'b');
ylabel('ISGs');
axis([0 tspan 1e-3 1e3])
subplot(6,1,4)
semilogy(Glu93Ala271_time,Glu93Ala271_Network_soln(:,4),'r');
hold on
semilogy(Glu93Gln271_time,Glu93Gln271_Network_soln(:,4),'g');
semilogy(Gln93Gln271_time,Gln93Gln271_Network_soln(:,4),'b');
ylabel('C3pro');
axis([0 tspan 1e-3 1e3])
subplot(6,1,5)
semilogy(Glu93Ala271_time,Glu93Ala271_Network_soln(:,5),'r');
hold on
semilogy(Glu93Gln271_time,Glu93Gln271_Network_soln(:,5),'g');
semilogy(Gln93Gln271_time,Gln93Gln271_Network_soln(:,5),'b');
ylabel('MAVSclv1');
axis([0 tspan 1e-3 1e3])
subplot(6,1,6)
semilogy(Glu93Ala271_time,Glu93Ala271_Network_soln(:,6),'r');
hold on
semilogy(Glu93Gln271_time,Glu93Gln271_Network_soln(:,6),'g');
semilogy(Gln93Gln271_time,Gln93Gln271_Network_soln(:,6),'b');
ylabel('MAVSclv2');
axis([0 tspan 1e-3 1e3])
xlabel('Time')
legend({'Glu^9^3Ala^2^7^1' 'Glu^9^3Gln^2^7^1' 'Gln^9^3Gln^2^7^1'}, ...
    'Position',[0.8205 0.3667 0.1634 0.1083]);

%ISG sensitivity analysis for protein activation rate
logfactor = 2;  %fine scaling factor for sensitivity analysis
k_act_vect = [k_act/logfactor^2 k_act/logfactor ...
    k_act k_act*logfactor k_act*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(k_act_vect)
    model_params_temp(1)=k_act_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    k_act_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))]; %time-integrate ISGs
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,1)
loglog(k_act_vect,k_act_sens(:,1),'ro-')
hold on
loglog(k_act_vect,k_act_sens(:,2),'go-')
loglog(k_act_vect,k_act_sens(:,3),'bo-')
xlabel('Protein activation rate')
ylabel('Time-integrated ISGs')
axis([min(k_act_vect) max(k_act_vect) 0.1 1000])

%ISG sensitivity analysis for MAVS polymerization rate
logfactor = 2;  %fine scaling factor for sensitivity analysis
k_poly_vect = [k_poly/logfactor^2 k_poly/logfactor ...
    k_poly k_poly*logfactor k_poly*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(k_poly_vect)
    model_params_temp(2)=k_poly_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    k_poly_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))]; %time-integrate ISGs
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,2)
loglog(k_poly_vect,k_poly_sens(:,1),'ro-')
hold on
loglog(k_poly_vect,k_poly_sens(:,2),'go-')
loglog(k_poly_vect,k_poly_sens(:,3),'bo-')
xlabel('MAVS polymerization rate')
ylabel('Time-integrated ISGs')
axis([min(k_poly_vect) max(k_poly_vect) 0.1 1000])

%ISG sensitivity analysis for protein degradation rate
logfactor = 2;  %fine scaling factor for sensitivity analysis
k_deg_vect = [k_deg/logfactor^2 k_deg/logfactor ...
    k_deg k_deg*logfactor k_deg*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(k_deg_vect)
    model_params_temp(3)=k_deg_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    k_deg_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))]; %time-integrate ISGs
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,3)
loglog(k_deg_vect,k_deg_sens(:,1),'ro-')
hold on
loglog(k_deg_vect,k_deg_sens(:,2),'go-')
loglog(k_deg_vect,k_deg_sens(:,3),'bo-')
xlabel('Protein degradation rate')
ylabel('Time-integrated ISGs')
axis([min(k_deg_vect) max(k_deg_vect) 0.1 1000])

%ISG sensitivity analysis for ISG transcription rate
logfactor = 2;  %fine scaling factor for sensitivity analysis
k_txn_vect = [k_txn/logfactor^2 k_txn/logfactor ...
    k_txn k_txn*logfactor k_txn*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(k_txn_vect)
    model_params_temp(4)=k_txn_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    k_txn_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))]; %time-integrate ISGs
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,4)
loglog(k_txn_vect,k_txn_sens(:,1),'ro-')
hold on
loglog(k_txn_vect,k_txn_sens(:,2),'go-')
loglog(k_txn_vect,k_txn_sens(:,3),'bo-')
xlabel('ISG transcription rate')
ylabel('Time-integrated ISGs')
axis([min(k_txn_vect) max(k_txn_vect) 0.1 1000])

%ISG sensitivity analysis for cleavage rate
logfactor = 2;  %finescaling factor for sensitivity analysis
k_clv_vect = [k_clv/logfactor^2 k_clv/logfactor ...
    k_clv k_clv*logfactor k_clv*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(k_clv_vect)
    model_params_temp(5)=k_clv_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    k_clv_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))];  %time-integrate ISGs
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,5)
loglog(k_clv_vect,k_clv_sens(:,1),'ro-')
hold on
loglog(k_clv_vect,k_clv_sens(:,2),'go-')
loglog(k_clv_vect,k_clv_sens(:,3),'bo-')
xlabel('Pro3C cleavage rate')
ylabel('Time-integrated ISGs')
axis([min(k_clv_vect) max(k_clv_vect) 0.1 1000])

%ISG sensitivity analysis for negative-feedback strength
logfactor = 10;  %coarse scaling factor for sensitivity analysis
phi_vect = [phi/logfactor^2 phi/logfactor ...
    phi phi*logfactor phi*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(phi_vect)
    model_params_temp(6)=phi_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    phi_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))];  %time-integrate ISGs
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,6)
loglog(phi_vect,phi_sens(:,1),'ro-')
hold on
loglog(phi_vect,phi_sens(:,2),'go-')
loglog(phi_vect,phi_sens(:,3),'bo-')
xlabel('Negative-feedback strength')
ylabel('Time-integrated ISGs')
axis([min(phi_vect) max(phi_vect) 0.1 1000])

%ISG sensitivity analysis for ratio of negative-feedback strength
logfactor = 10;  %scaling factor for sensitivity analysis
gamma_vect = [gamma/logfactor^2 gamma/logfactor ...
    gamma gamma*logfactor gamma*logfactor^2];
model_params_temp = Model_Params;
for i=1:length(gamma_vect)
    model_params_temp(7)=gamma_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    gamma_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))];
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,7)
loglog(gamma_vect,gamma_sens(:,1),'ro-')
hold on
loglog(gamma_vect,gamma_sens(:,2),'go-')
loglog(gamma_vect,gamma_sens(:,3),'bo-')
xlabel('Ratio of negative-feedback strengths')
ylabel('Time-integrated ISGs')
axis([min(gamma_vect) max(gamma_vect) 0.1 1000])

%ISG sensitivity analysis for filament size
filament_vect = linspace(400,1000,28);
model_params_temp = Model_Params;
for i=1:length(filament_vect)
    model_params_temp(8)=filament_vect(i);
    [Glu93Ala271_time_temp,Glu93Ala271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Ala271');
    [Glu93Gln271_time_temp,Glu93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Glu93Gln271');
    [Gln93Gln271_time_temp,Gln93Gln271_Network_soln_temp] = ode15s(@NetworkODE, ...
        [0 tspan],Network_IC,[],model_params_temp,'Gln93Gln271');
    filament_sens(i,:)=[trapz(Glu93Ala271_time_temp, ...
        Glu93Ala271_Network_soln_temp(:,3)) trapz(Glu93Gln271_time_temp, ...
        Glu93Gln271_Network_soln_temp(:,3)) trapz(Gln93Gln271_time_temp, ...
        Gln93Gln271_Network_soln_temp(:,3))];
    clear Glu93Ala271_time_temp Glu93Ala271_Network_soln_temp Glu93Gln271_time_temp ...
        Glu93Gln271_Network_soln_temp Gln93Gln271_time_temp Gln93Gln271_Network_soln_temp
end
figure(2)
subplot(3,3,8)
semilogy(filament_vect,filament_sens(:,1),'ro-')
hold on
semilogy(filament_vect,filament_sens(:,2),'go-')
semilogy(filament_vect,filament_sens(:,3),'bo-')
xlabel('MAVS filament size')
ylabel('Time-integrated ISGs')
axis([400 1000 0.1 1000])
legend({'Glu^9^3Ala^2^7^1' 'Glu^9^3Gln^2^7^1' 'Gln^9^3Gln^2^7^1'}, ...
    'Position',[0.7093 0.1661 0.1634 0.1083]);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dNetwork_dt = NetworkODE(t,Network,params,model)

param_cell = num2cell(params);
[k_a,k_p,k_d,k_t,k_c,fb,g,n] = param_cell{:};

Network_cell = num2cell(Network);
[MAVS,polyMAVS,ISGs,C3pro,MAVSclv1,MAVSclv2] = Network_cell{:};


if strcmpi(model,'Glu93Ala271')
    dMAVS_dt = k_a - k_p*MAVS.^n;
    dpolyMAVS_dt = k_p*MAVS.^n - k_d/n*polyMAVS;
    dISGs_dt = k_t/n*polyMAVS - k_d*ISGs;
    dC3pro_dt = k_a;
    dMAVSclv1_dt = 0;   %Glu93Gln271 cleavage site is mutated
    dMAVSclv2_dt = 0;   %polymorphism is absent
end

if strcmpi(model,'Glu93Gln271')
    dMAVS_dt = k_a - k_p*(fb/(fb+MAVSclv1))*MAVS.^n - k_c*C3pro*MAVS;
    dpolyMAVS_dt = k_p*(fb/(fb+MAVSclv1))*MAVS.^n - k_d/n*polyMAVS;
    dISGs_dt = k_t/n*polyMAVS - k_d*ISGs;
    dC3pro_dt = k_a;
    dMAVSclv1_dt = k_c*C3pro*MAVS;
    dMAVSclv2_dt = 0;   %polymorphism is absent
end

if strcmpi(model,'Gln93Gln271')
    dMAVS_dt = k_a - k_p*(fb/(fb+MAVSclv2+MAVSclv1/g))*MAVS.^n - ...
        2*k_c*C3pro*MAVS;
    dpolyMAVS_dt = k_p*(fb/(fb+MAVSclv2+MAVSclv1/g))*MAVS.^n - ...
        k_d/n*polyMAVS;
    dISGs_dt = k_t/n*polyMAVS - k_d*ISGs;
    dC3pro_dt = k_a;
    dMAVSclv1_dt = k_c*C3pro*MAVS;
    dMAVSclv2_dt = k_c*C3pro*MAVS;
end

dNetwork_dt = [dMAVS_dt; dpolyMAVS_dt; dISGs_dt; dC3pro_dt; dMAVSclv1_dt; ...
    dMAVSclv2_dt];

return