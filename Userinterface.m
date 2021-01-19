close all
fig = uifigure;
trans_lbl = uilabel(fig,'Position',[210 380 120 22],'Text','TRANSMITTER', 'FontSize',16);
cp_lbl = uilabel(fig,'Position',[218 340 120 22],'Text','Enter CP length', 'FontSize',14);
cp_txt = uitextarea(fig,'Position',[220 320 100 22],'Value',' ');
blk_dd = uidropdown(fig,'Position',[160 280 230 30],'Items',{'Select number of blocks for GFDM','2','4','8','16','32'},...
     'Value','Select number of blocks for GFDM');
ts_btn = uibutton(fig,'Position',[200 240 150 30],'Text','Plot Spectrum','ButtonPushedFcn', @(ts_btn,event) plotSpectrum(ts_btn,str2num(cell2mat(cp_txt.Value)),256/(str2double(get(blk_dd, 'Value')))));
rcv_lbl = uilabel(fig,'Position',[230 190 120 22],'Text','RECEIVER', 'FontSize',16);
blk_lbl = uilabel(fig,'Position',[200 160 250 22],'Text','Drag Slider to Select SNR', 'FontSize',14);
snr_slider = uislider(fig,'Position',[150 150 250 30],'Limits',[-10 10]);
rcv_btn = uibutton(fig,'Position',[175 80 200 30],'Text','See Receiver Performance','ButtonPushedFcn',@(rcv_btn,event) Rec_performance(rcv_btn, str2num(cell2mat(cp_txt.Value)),256/(str2double(get(blk_dd, 'Value'))),ceil(snr_slider.Value)));

function plotSpectrum(ts_btn, cp_, sym_num)
cp = cp_; no_of_sym = sym_num;
trans = transmitter();
trans.mapper();
ofdm_sig = trans.modulate('OFDM',cp,1);
gfdm_sig = trans.modulate('GFDM',cp,no_of_sym);
close all
figure 
f1 = linspace(-trans.no_of_total_subcarriers/2, trans.no_of_total_subcarriers/2, 2*length(ofdm_sig)+1);
f1 = f1(1:end-1)'; f1 = normalize(f1,'range', [-1 1]);
f2 = linspace(-trans.no_of_total_subcarriers/2, trans.no_of_total_subcarriers/2, 2*length(gfdm_sig)+1); 
f2 = f2(1:end-1)'; f2 = normalize(f2, 'range', [-1 1]);
plot(f1, mag2db(flip(abs(fft(ofdm_sig, 2*length(ofdm_sig)))))/2, 'b');
title('OFDM Spectrum'), 
xlabel('normalize frequency'), ylabel('Power Spectral Density [dB]');ylim([-10 30])
figure
plot(f2, mag2db(flip(abs(fft(gfdm_sig, 2*length(gfdm_sig)))))/2, 'b');
title('GFDM Spectrum'),
xlabel('normalize frequency'), ylabel('Power Spectral Density [dB]'); ylim([-10 30])
end

function Rec_performance(rcv_btn, cp_, sym_num, snr)
cp = cp_; no_of_sym = sym_num; user_snr = snr;
close all;
[oser, ober] = OFDMreceiver(user_snr,cp);
[gser, gber] = GFDMreceiver(user_snr,cp,no_of_sym);
fprintf('OFDM BER = %f at SNR = %d \n',ober,snr); fprintf('OFDM SER = %f at SNR = %d\n',oser,snr);
fprintf('GFDM BER = %f at SNR = %d \n',gber,snr); fprintf('GFDM SER = %f at SNR = %d\n',gser, snr);
end