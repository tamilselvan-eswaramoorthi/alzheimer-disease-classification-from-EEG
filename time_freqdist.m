function [Cxx,AF,wind] = time_freqdist(in_signal,sig_len,vec,tau,lambda,ansig)




if nargin<6
    ansig = 1;
end
if ansig
    in_signal = hilbert(real(in_signal));
end

% Wigner-Ville distribution
wvx=tfrwv(in_signal(:),1:length(in_signal),sig_len);
% Ambiguity function
AF=fft(ifft(wvx).');clear wvx

if rem(size(AF,2),512)~=0 & rem(size(AF,2),1024)~=0 & rem(size(AF,2),4096)~=0  & rem(size(AF,2),256)~=0
    AF=AF.';
end
% k = 1:size(AF,2);

[Crow,Ccol]=size(AF);
dy=[-Ccol/2:Ccol/2-1]/(Ccol/2);
if Crow/2- fix(Crow/2)==0
    dx=[-Crow/2:Crow/2-1]/(Crow/2);
else
    dx=[-Crow/2-1/2:Crow/2-3/2]/(Crow/2);
end
[in_signal,y]=meshgrid(dy,dx);
tau1=in_signal/tau;
v1=y/vec;

beta = 1; gamma = 1; alfa = 0; r = 0; 
mu=  (tau1.^2.*(v1.^2).^alfa + (tau1.^2).^alfa .* v1.^2 +2.*r*((tau1.*v1).^beta).^gamma);
clear tau1 v1 x y
wind = exp(-pi*( mu.^2).^lambda);
wind = wind/max(max(wind));
Asmooth = (AF.*fftshift(wind)).';
Cxx = real(ifft(fft(Asmooth).').');