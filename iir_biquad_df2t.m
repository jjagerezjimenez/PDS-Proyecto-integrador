%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCION: iir_biquad_df2t
%  Procesa 1 muestra con una celda de 2do orden en realización directa 
%  transpuesta (DF2T).
%  
%  Entrada:
%   x_in   : muestra de entrada
%   sosVec : [b0 b1 b2 a0 a1 a2] coeficientes de la sección
%   states : [w1 w2] registros de demora
%
%  Salida:
%   y_out  : muestra de salida
%   states : estados actualizados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_out, states] = iir_biquad_df2t(x_in, sosVec, states)
    % Desempaquetamos los coeficientes
    b0 = sosVec(1); b1 = sosVec(2); b2 = sosVec(3);
    a0 = sosVec(4); a1 = sosVec(5); a2 = sosVec(6);

    % Estados anteriores w1 y w2
    w1 = states(1);
    w2 = states(2);

    % Operación DF2T:
    % w = x_in - a1*w1 - a2*w2
    % y = b0*w + b1*w1 + b2*w2
    w  = x_in - a1*w1 - a2*w2;
    y  = b0*w + b1*w1 + b2*w2;

    % Actualizamos los estados: w2 <- w1, w1 <- w
    states(2) = w1; 
    states(1) = w;

    % Si a0 != 1, normalizamos la salida
    if a0~=1
        y = y/a0;
    end

    % Devolvemos la salida y los estados
    y_out = y;
end
