%% FYS-6607 Thin-films and waveguides-course, 2021
% Tapio Niemi and Janne Simonen, 
% tapio.niemi@tuni.fi

%%   Calculates fields for a multi-layer stack

% Formulation follows closely the approach by Yeh and Yariv
% Photonics, ISBN 13: 9780195179460, Oxford University Press
% or
% Optical Waves in Layered Media, ISBN-13: 978-0471731924
% Wiley Series in Pure and Applied Optics

%% Parameters
% Input
% n1 = refractive index of incidence side material
% n2 = refractive index after the thin film stack
% n = refractive indices of the layers of the stack, vector
% L = thicknesses of the layers of the stack, vector [in nanometers]
% lambda = wavelength [in nanometers]
% theta = incidence angle [in degrees]
% polarization = 'TE' or 'TM'

% Output
% z   = length along the stack
% nz  = refractive index profile along the stack
% Et   = Tangential electric field 
% En   = Normal electric field (=0 for TE)
% Stot = Poynting vector


function [z, nz, E0, E_t, E_n, Stot] = tmatrix_fields(n1,n2,n,L,lambda,theta,polarization)
        
% First calculate the reflection of the stack
% This is needed as an input
[refl,~,~,~]=tmatrix_filmstack(n1,n2,n,L,lambda,theta,polarization);

% Refractive index vector including surrounding materials
% Thickness vector including the sub-/superstrate
n=[n1 n n2];
cthickness=2*lambda;
L=[-cthickness L cthickness]; %+-cthickness into the sub-/superstrate
N=length(n);

% Resolution of z-length
res=lambda/2^6;

% wavenumbers
kt=2*pi./lambda*n(1)*sin(theta*pi/180);
kn=sqrt((2*pi./lambda.*n).^2-kt.^2);
k=2*pi./lambda.*n;

% Note the units!! 
% Should be in [m]
omega=2*pi*299792458/lambda;
u0=4*pi*1e-7;

% Ensure correct sign of the imaginary part
% Especially in a case of total internal reflection
kn=real(kn)-1j*abs(imag(kn))

% This part depends on the polarization
% Reflectance and transmittance coefficients of each interface,
% vectors of length N-1
    if strcmp(polarization,'TE')
        
        r=(kn(1:N-1)-kn(2:N))./...
        (kn(1:N-1)+kn(2:N));
        t=(2*kn(1:N-1))./...
            (kn(1:N-1)+kn(2:N));

            elseif strcmp(polarization,'TM')
                
        r=(n(1:N-1).^2.*kn(2:N)-n(2:N).^2.*kn(1:N-1))./...
            (n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));
        t=(2*n(2:N).*n(1:N-1).*kn(1:N-1))./...
            (n(2:N).^2.*kn(1:N-1)+n(1:N-1).^2.*kn(2:N));

    else
        error('Polarization must be TE or TM!')
    end
    
% Input amplitude E0=1 and reflected is r*E0 
% Field amplitude at the first interface is then
A=[1;refl];

% Generate z-axis from -L to 0
Nz=abs(floor(L(1)/res));
z=linspace(L(1),0,Nz);

% E-field is a sum of forward and backward waves in each layer
% E-field from -L to 0
E0=A(1,1)*exp(-1j*kn(1)*z)+A(2,1)*exp(1j*kn(1)*z);
    
    if strcmp(polarization,'TE')
        Ecladding1_t=E0;
        Ecladding1_n=zeros(size(Ecladding1_t));
        
    else % TM polarization
            
        Ecladding1_n=kt/k(1)*E0;
        Ecladding1_t=kn(1)/k(1)*E0;
    
    end
        
    E_t=Ecladding1_t;
    E_n=Ecladding1_n;
             
    % Poynting vector from -L to 0. 
    S_A=1/(2*omega*u0)*real(kn(1))*abs(A(1,1))^2*exp(2*imag(kn(1))*z);
    S_B=1/(2*omega*u0)*real(kn(1))*abs(A(2,1))^2*exp(-2*imag(kn(1))*z);
    S_AB=1/(omega*u0)*imag(kn(1))*imag(A(1,1)*conj(A(2,1))*exp(-2j*real(kn(1))*z));
    Stot=S_A-S_B-S_AB;
    
    Stot_A=S_A; % Forward flow
    Stot_B=S_B; % Backward flow
    Stot_AB=S_AB; % Coupling due to loss
    
    % Generate refractive index vector from -L to 0
    nz=n(1)*ones(1,length(z));
                        
% Loop the other layers from 2 to final
    
    for M=2:length(L)
    
        % Transmission matrix through the interface M
        T=1/t(M-1)*[1 r(M-1);r(M-1) 1];
        
        % Field amplitudes after the interface M
       % Tinv=1/(1-r(M-1)^2)/t(M-1)*[1 -r(M-1);-r(M-1) 1]
        A=inv(T)*A;
        %A=Tinv*A;
        
        % Points for the layer M from 0 to L
        Nz=abs(floor(L(M)/res));
        z_LayerM=linspace(res,L(M),Nz);
        
        %Update z-vector
        z=[z z(end)+z_LayerM];
        A
        % Field for the layer M from 0 to L
        E0_M=A(1,1)*exp(-1j*kn(M)*z_LayerM)+A(2,1)*exp(1j*kn(M)*z_LayerM);
 
        if strcmp(polarization,'TE')
            
            E_tM=E0_M;
            E_nM=zeros(size(E_tM));
        
                else % TM polarization
            
            E_nM=kt/k(M)*E0_M;
            E_tM=kn(M)/k(M)*E0_M;
    
        end
        %Update E-field
        E0=[E0 E0_M];
        E_t=[E_t E_tM];
        E_n=[E_n E_nM];
                
        % Poynting vector from 0 to L in layer M
        S_A=1/(2*omega*u0)*real(kn(M))*abs(A(1,1))^2*exp(2*imag(kn(M))*z_LayerM);
        S_B=1/(2*omega*u0)*real(kn(M))*abs(A(2,1))^2*exp(-2*imag(kn(M))*z_LayerM);
        S_AB=1/(2*omega*u0)*2*imag(kn(M))*imag(A(1,1)*conj(A(2,1))*exp(-2j*real(kn(M))*z_LayerM));
            
        S_layerM=S_A-S_B-S_AB;
                     
        Stot_A=[Stot_A S_A];
        Stot_B=[Stot_B S_B];
        Stot_AB=[Stot_AB S_AB];
        Stot=[Stot S_layerM];
        
        % Field amplitudes after the propagation
        Pinv=[exp(-j*kn(M)*L(M)) 0 ; 0 exp(j*kn(M)*L(M))];
        A=Pinv*A; %E0(:,end);
                
        % Update refractive index vector
        nz=[nz n(M)*ones(1,length(z_LayerM))];
    end  
        % Total power normalized to incident power
        Stot=Stot/Stot_A(1);
        Stot_A(1)
end