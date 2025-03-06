%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCION AUXILIAR: fftPlot
%   Realiza la FFT de "sig" y grafica la magnitud en dB desde 0..freqMax
%   para señal real.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fftPlot(sig, fs, freqMax)
    % Calcula el tamaño de FFT como potencia de 2
    Nfft = 2^nextpow2(length(sig));   
    % FFT de la señal con zero-padding
    SIG = fft(sig, Nfft);
    % Magnitud
    magSIG = abs(SIG);
    % Magnitud en dB, evitando log(0)
    magSIGdB = 20*log10(magSIG + 1e-12);

    % Eje de frecuencias
    f = (0:Nfft-1)*(fs/Nfft);
    halfRange = 1:(Nfft/2); 
    fPlot = f(halfRange);
    magPlot = magSIGdB(halfRange);

    % Recortamos a freqMax si es menor que fs/2
    idxMax = find(fPlot <= freqMax, 1, 'last');
    if isempty(idxMax)
        idxMax = length(fPlot);
    end

    % Graficamos
    plot(fPlot(1:idxMax)/1000, magPlot(1:idxMax),'LineWidth',1);
    xlabel('Frecuencia [kHz]'); ylabel('Magnitud [dB]');
    grid on;
end