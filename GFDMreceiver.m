trans = transmitter();
trans.mapper('QPSK');
trans.modulate('GFDM',8,64);
GFDM_signal = trans.transmitted_signal;
sen = trans.binary_source;
g = trans.filter_pulse;
snr_dB = -10:9;
errors   = zeros(length(snr_dB), 1);
s_ = zeros(trans.no_of_total_subcarriers,trans.no_of_subsymbols,trans.no_of_blocks);
sf = zeros(trans.no_of_subcarriers,trans.no_of_subsymbols,trans.no_of_blocks);
s2 = zeros(trans.no_of_subcarriers*trans.no_of_subsymbols,trans.no_of_blocks);
s3 = zeros(1,trans.no_of_subcarriers*trans.no_of_subsymbols*trans.no_of_blocks);
n = 1/sqrt(2)*(randn(1,length(GFDM_signal))+ 1i*randn(1,length(GFDM_signal)));
for si = 1:length(snr_dB)
    rcv_signal = GFDM_signal + 10^(-snr_dB(si)/20)*n; %pass through channel
    x = circshift(rcv_signal'.*exp(1i*2*pi*0*(1:length(rcv_signal))'), 0);
    y = reshape(x,[],trans.no_of_blocks);
    y = y(trans.len_of_CP+1:end,:); % remove CP
    IQ = zeros(trans.no_of_total_subcarriers,trans.no_of_subsymbols,trans.no_of_blocks);
    D = zeros(trans.no_of_total_subcarriers,trans.no_of_subsymbols,trans.no_of_blocks);
    s = zeros(trans.no_of_total_subcarriers*trans.no_of_subsymbols,trans.no_of_blocks);
    for b = 1:trans.no_of_blocks
        yhat = y(:,b);
        G = conj(fft(g));
        L = length(G)/trans.no_of_subsymbols;
        xhat = fft(yhat);
        Dhat = zeros(trans.no_of_total_subcarriers,trans.no_of_subsymbols);
       for k = 1:trans.no_of_total_subcarriers
           carrier = circshift(xhat,ceil(L*trans.no_of_subsymbols/2) - trans.no_of_subsymbols*(k-1));
           carrier = fftshift(carrier(1:L*trans.no_of_subsymbols));
           carrierMatched = carrier .* G;
           dhat = ifft(sum(reshape(carrierMatched,trans.no_of_subsymbols,L),2)/L);
           D(k,:,b) = dhat;
       end
         s(:,b) = reshape(D(:,:,b), numel(D(:,:,b)),1);
         s_(:,:,b) = reshape(s(:,b),trans.no_of_total_subcarriers,trans.no_of_subsymbols);
         sf(:,:,b) = s_(51:trans.no_of_subcarriers+50,:,b); %remove null carriers
         s2(:,b) = reshape(sf(:,:,b),[],1);
    end
    s3 = reshape(s2,1,[]);
    sk = qamdemod(s3, 2^2,'gray','OutputType','bit','UnitAveragePower',true);
    sx = reshape(sk,1,[]);
    errors(si) = length(find(sx' - sen'));
    
end
ber = errors/length(sx);
plot(snr_dB, log10(ber), '-x');
snrlin=10.^(snr_dB./10);
formula_ber = 0.5*erfc(sqrt(snrlin));
hold on
plot(snr_dB, log10(formula_ber),'-o');
