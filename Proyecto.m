%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROYECTO FDM 4 CANALES, IIR ORDEN 8 (CELDAS DF2T) + GRÁFICAS FFT
%
% Se generan gráficas de:
%  1) FFT de cada señal original (0..5 kHz)
%  2) FFT de cada señal tras interpolación
%  3) FFT de la señal en banda ancha
%  4) FFT de los filtros FIR
%  5) FFT de las señales filtradas (antes de la suma)
%  6) FFT de las señales recuperadas
%
% Y se demultiplexa con un IIR de orden 8 -> 4 celdas de 2º orden en cascada,
% en realización directa transpuesta, procesando muestra a muestra.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;  % Limpia todas las variables del workspace
close all;  % Cierra todas las ventanas de figuras abiertas
clc;    % Limpia la ventana de comandos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTE 1: Lectura de audios y parámetros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hector, fs1]  = audioread('hector8k.wav');  % Lee el archivo hector8k.wav
[santi,  fs2]  = audioread('santi8k.wav');   % Lee el archivo santi8k.wav
[juanjo, fs3]  = audioread('juanjo8k.wav');  % Lee el archivo juanjo8k.wav
[ernesto,fs4]  = audioread('ernesto8k.wav'); % Lee el archivo ernesto8k.wav

% Verifica que todos los archivos estén a 8 kHz
if any([fs1, fs2, fs3, fs4] ~= 8000)
    error('Algún archivo no está a 8 kHz');
end

% Asegura que sean vectores fila
hector  = hector(:).';  
santi   = santi(:).';
juanjo  = juanjo(:).';
ernesto = ernesto(:).';

% Calcula la longitud de cada archivo
N1 = length(hector); 
N2 = length(santi);
N3 = length(juanjo);
N4 = length(ernesto);

% Encuentra la longitud mínima para recortar todos igual
Nmin = min([N1,N2,N3,N4]);

% Recorta cada archivo a la misma longitud Nmin
hector  = hector(1:Nmin);
santi   = santi(1:Nmin);
juanjo  = juanjo(1:Nmin);
ernesto = ernesto(1:Nmin);

% Almacena las 4 señales en la matriz "voz", filas = canales
voz = [hector; santi; juanjo; ernesto];

nCanales = 4;       % Número de canales
fs_bb   = 8000;     % Frecuencia de muestreo banda base
fs_FDM  = 328000;   % Frecuencia de muestreo banda ancha
M       = fs_FDM/fs_bb;  % Factor 41

fprintf('Trabajando con %d canales, cada uno con %d muestras a 8 kHz.\n',...
    nCanales, Nmin);

%% 1) FFT de cada señal ORIGINAL (0..5 kHz)
figure('Name','FFT de cada señal de audio original');  % Crea una figura
for c = 1:nCanales
    subplot(nCanales,1,c);  % Subplot para cada canal
    fftPlot(voz(c,:), fs_bb, 6000);  % Llama la función fftPlot hasta 6 kHz
    title(sprintf('Canal %d original - FFT (0..5kHz)', c)); % Título
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTE 2: Filtros FIR pasa banda, Overlap&Add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definimos las bandas de cada canal
iniciosBanda = [60000, 64000, 68000, 72000];
finesBanda   = [64000, 68000, 72000, 76000];

L    = 1024;   % Longitud del FIR
Nblk = 4096;   % Tamaño del bloque Overlap&Add
Hf_cell = cell(1,nCanales);  % Guardará la FFT de cada FIR

% Para graficar la FFT de los filtros, guardamos los coeficientes
coefFiltro = cell(1, nCanales);

% Bucle para diseñar y almacenar la FFT de cada FIR
for c = 1:nCanales
    f1 = iniciosBanda(c);
    f2 = finesBanda(c);
    fc1 = f1/(fs_FDM/2);
    fc2 = f2/(fs_FDM/2);

    % Filtro FIR pasa banda, ventana Kaiser(4)
    htemp = fir1(L-1, [fc1, fc2], 'bandpass', kaiser(L,4));

    % Normalizar la ganancia en la frecuencia central
    fcMid = (f1+f2)/2;
    [HtempFreq, Wfreq] = freqz(htemp, 1, 2048, fs_FDM);
    [~, idxMid] = min(abs(Wfreq - fcMid));
    gainMid = abs(HtempFreq(idxMid));
    if gainMid~=0
        htemp = htemp / gainMid;
    end

    % Calcula la FFT del filtro con zero-padding hasta Nblk
    Hf_cell{c} = fft(htemp, Nblk);

    % Guarda los coef para graficar su respuesta en frecuencia
    coefFiltro{c} = htemp;
end

%% 4) GRÁFICAS FFT DE LOS FILTROS
figure('Name','FFT de los filtros FIR (4 kHz)');
for c = 1:nCanales
    subplot(nCanales,1,c);
    HH = fft(coefFiltro{c}, 4096);         % FFT de los coeficientes
    ejeFrecuencia = (0:4095)*(fs_FDM/4096);% Eje de frecuencia
    magH = 20*log10(abs(HH)+1e-12);        % Magnitud en dB
    plot(ejeFrecuencia, magH, 'LineWidth',1);
    xlim([0, 120e3]);  % Limitamos hasta 120 kHz
    grid on;
    xlabel('Frecuencia [Hz]'); ylabel('dB');
    title(sprintf('Filtro canal %d (FIR) - FFT mag', c));
end

% Preparar buffers para Overlap&Add
tempBuff  = zeros(nCanales, Nblk);    % Buffer de tamaño Nblk por canal
overlap   = zeros(nCanales, L-1);     % Para manejar el solape (L-1)
posBuff   = zeros(1, nCanales);       % Índice dentro del buffer para cada canal

% Matriz para las señales filtradas (antes de sumarlas)
canalesFiltrados = cell(1, nCanales);
for c = 1:nCanales
    canalesFiltrados{c} = zeros(1, Nmin*M);
end

% Para graficar las señales luego de la interpolación de ceros
senialesInterpoladas = cell(1,nCanales);
for c = 1:nCanales
    senialesInterpoladas{c} = zeros(1, Nmin*M);
end

indiceFDM = 1;  % Índice global en la señal banda ancha

% Bucle principal para cada muestra en banda base
for n = 1:Nmin
    for c = 1:nCanales
        x_c = voz(c,n);         % Muestra actual del canal c
        factorSigno = 1;
        if mod(c,2)==1          % Canales impares => alternan signo
            factorSigno = (-1)^(n-1);
        end

        % Interpola x41 => 1 muestra real y 40 ceros
        for k=1:M
            if k==1
                muestraEntrada = x_c*factorSigno;
            else
                muestraEntrada = 0;
            end

            % Guardamos la señal interpolada para graficar su FFT
            idxInterp = (n-1)*M + k;
            senialesInterpoladas{c}(idxInterp) = muestraEntrada;

            % Llenamos el buffer para Overlap&Add
            posBuff(c) = posBuff(c) + 1;
            tempBuff(c, posBuff(c)) = muestraEntrada;

            % Si llenamos el buffer (Nblk), procesamos Overlap&Add
            if posBuff(c) == Nblk
                Xf = fft(tempBuff(c,:), Nblk);   % FFT del bloque
                Yf = Xf .* Hf_cell{c};          % Multiplica en freq
                y  = ifft(Yf, Nblk);            % IFFT => salida

                % Solape
                y(1:(L-1)) = y(1:(L-1)) + overlap(c,:);
                overlap(c,:) = y(Nblk-(L-2):Nblk);

                % Ubicamos el resultado en canalesFiltrados
                idxStart = indiceFDM - (Nblk-1);
                if idxStart < 1
                    idxStart = 1;
                end
                idxEnd = indiceFDM;
                segLen = idxEnd - idxStart + 1;
                if idxEnd > (Nmin*M)
                    idxEnd = Nmin*M;
                    segLen = idxEnd - idxStart + 1;
                end

                canalesFiltrados{c}(idxStart:idxEnd) = ...
                   canalesFiltrados{c}(idxStart:idxEnd) + y(1:segLen);

                % Reiniciamos buffer
                tempBuff(c,:) = 0;
                posBuff(c)    = 0;
            end
        end
    end
    % Avanzamos indiceFDM en M muestras
    indiceFDM = indiceFDM + M;
end

% Sumamos las 4 salidas filtradas => señal banda ancha FDM
fdm_signal = zeros(1, Nmin*M);
for c = 1:nCanales
    fdm_signal = fdm_signal + canalesFiltrados{c};
end

fprintf('Señal FDM generada: %d muestras a Fs=%d\n',...
    length(fdm_signal), fs_FDM);

%% 2) FFT SEÑALES ORIGINALES LUEGO DE INTERPOLAR
figure('Name','FFT tras Interpolar (canales)');
for c = 1:nCanales
    subplot(nCanales,1,c);
    fftPlot(senialesInterpoladas{c}, fs_FDM, 160000);
    title(sprintf('Canal %d - FFT tras Interpolar', c));
end

%% 3) FFT SEÑAL EN BANDA ANCHA
figure('Name','Señal FDM banda ancha');
fftPlot(fdm_signal, fs_FDM, 120000);
title('Señal FDM banda ancha');

%% 5) FFT DE LAS SEÑALES FILTRADAS (antes de sumar)
figure('Name','FFT de las señales filtradas (Overlap&Add)');
for c = 1:nCanales
    subplot(nCanales,1,c);
    fftPlot(canalesFiltrados{c}, fs_FDM, 90000);
    title(sprintf('Canal %d filtrado (Overlap&Add)', c));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTE 3: DEMULTIPLEXADO - IIR ORDEN 8 (4 celdas DF2T) MUESTRA A MUESTRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ordenBP = 8;  % => 4 celdas de 2º orden en cascada
sos_cell = cell(1,nCanales);

% Diseño de filtros IIR (Butterworth orden 8)
for c = 1:nCanales
    f1 = iniciosBanda(c);
    f2 = finesBanda(c);
    fc1 = f1/(fs_FDM/2);
    fc2 = f2/(fs_FDM/2);

    [b_iir, a_iir] = butter(ordenBP, [fc1 fc2]);

    % Normalizar en frecuencia central
    fcMid = (f1+f2)/2;
    [Hiir, Wiir] = freqz(b_iir, a_iir, 2048, fs_FDM);
    [~, idxMidIIR] = min(abs(Wiir - fcMid));
    gMid = abs(Hiir(idxMidIIR));
    if gMid~=0
        b_iir = b_iir / gMid;
    end

    sos_cell{c} = tf2sos(b_iir, a_iir);
end

% Estructura con sos y estados
dfiltIIR = cell(1,nCanales);
for c = 1:nCanales
    numSec = size(sos_cell{c},1); % para orden 8 => 4 secciones
    dfiltIIR{c}.sos = sos_cell{c};
    dfiltIIR{c}.estados = zeros(numSec,2); % cada sección 2 estados
end

% Procesamos muestra a muestra
senialesDemux = zeros(nCanales, length(fdm_signal));

for n = 1:length(fdm_signal)
    muestraEntrada = fdm_signal(n);
    for c = 1:nCanales
        ytemp = muestraEntrada;
        for sSec = 1:size(dfiltIIR{c}.sos,1)
            % Llamamos la función iir_biquad_df2t
            [ytemp, dfiltIIR{c}.estados(sSec,:)] = iir_biquad_df2t(...
                ytemp, ...
                dfiltIIR{c}.sos(sSec,:), ...
                dfiltIIR{c}.estados(sSec,:) ...
                );
        end
        senialesDemux(c,n) = ytemp;
    end
end

% Deshacer desplazamiento en impares
for c = 1:nCanales
    if mod(c,2)==1
        for n = 1:length(senialesDemux(c,:))
            senialesDemux(c,n) = senialesDemux(c,n)*((-1)^(n-1));
        end
    end
end

% Decimación /41 => fs=8kHz
recuperadas = cell(1,nCanales);
for c = 1:nCanales
    recuperadas{c} = senialesDemux(c, 1:M:end);
end

%% 6) FFT DE LAS SEÑALES RECUPERADAS
figure('Name','FFT de las señales recuperadas');
for c = 1:nCanales
    subplot(nCanales,1,c);
    fftPlot(recuperadas{c}, fs_bb, 6000);
    title(sprintf('Canal %d recuperado - FFT (0..5kHz)', c));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTE 4: Comparación, ganancia y reproducción
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gananciaRec = 41.0;  % factor de ganancia en recuperados

% Graficamos original vs recuperada
t_bb = (0:Nmin-1)/fs_bb;
figure('Name','Comparación: Original vs Recuperada (orden IIR=8)');
for c = 1:nCanales
    subplot(nCanales,1,c);
    RecSig = recuperadas{c} * gananciaRec;  % Aplicamos ganancia
    Lr = length(RecSig);
    Lm = min(Nmin, Lr);

    plot(t_bb(1:Lm), voz(c,1:Lm), 'b'); hold on;
    plot(t_bb(1:Lm), RecSig(1:Lm), 'r--');
    xlabel('Tiempo [s]'); ylabel('Amplitud');
    title(sprintf('Canal %d: Original vs Recuperada (Ganancia=%.1f, IIR=8)', ...
        c, gananciaRec));
    legend('Original','Recuperada');
    grid on;
end

%% 7) FFT DE LAS SEÑALES RECUPERADAS CON GANANCIA
figure('Name','FFT de las señales recuperadas con ganancia');
for c = 1:nCanales
    subplot(nCanales,1,c);
    RecSig = recuperadas{c} * gananciaRec;
    fftPlot(RecSig, fs_bb, 10000);
    title(sprintf('Canal %d recuperado con ganancia - FFT (0..5kHz)', c));
end

% Reproducción final
disp('Reproduciendo audios...');
for c = 1:nCanales
    durOri = length(voz(c,:))/fs_bb;
    sound(voz(c,:), fs_bb);          % Reproduce original
    pause(durOri + 0.5);

    senialSalida = recuperadas{c} * gananciaRec; % Recuperada c/ganancia
    durRec = length(senialSalida)/fs_bb;
    sound(senialSalida, fs_bb);
    pause(durRec + 0.5);
end

disp('Fin del proceso FDM.');






