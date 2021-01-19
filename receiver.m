%function receiver()
trans = transmitter();
trans.mapper('QPSK');
trans.modulate('OFDM',0,1);
OFDM_signal = trans.transmitted_signal;
%trans.modulate('GFDM',0,128);
%GFDM_signal = trans.transmitted_signal;
snr_dB = [-10:9];
errors = zeros(1,length(snr_dB));
berD = zeros(1,length(snr_dB));
n = 1/sqrt(2)*(randn(1,length(OFDM_signal))+ 1i*randn(1,length(OFDM_signal))); %Guassian noise

for i = 1:length(snr_dB)
    r_OFDM_signal = OFDM_signal + 10^(-snr_dB(i)/20)*n;
    r_OFDM_signal_matrix = reshape(r_OFDM_signal,trans.no_of_total_subcarriers,[]);
    r_OFDM_signal_ = r_OFDM_signal_matrix(:,trans.len_of_CP+1:end); %remove CP
    fft_OFDM_signal = fft(r_OFDM_signal_);
    fft_OFDM_signal_ = fft_OFDM_signal(51:trans.no_of_subcarriers+50,:);   % remove null subcarriers
    rcv_maped_signal = reshape(fft_OFDM_signal_,1,[]);
    rcv_bin_data_ = qamdemod(rcv_maped_signal,2^trans.no_of_bits_per_symbol,'OutputType','bit','UnitAveragePower',true);
    rcv_bin_data = reshape(rcv_bin_data_,1,[]);
    errors(i) = length(find(trans.binary_source-rcv_bin_data));
    %berD(i) = berawgn(snr_dB(i),'qam',16);
end
ber = errors/trans.no_of_bin_data;
ser = ber * trans.no_of_bits_per_symbol;
snrlin=10.^(snr_dB./10);
formula_ber = 0.5 * erfc(sqrt(snrlin));
plot(snr_dB, log10(ber), '-x');
hold on
plot(snr_dB, log10(formula_ber),'-o');

