classdef transmitter < handle
    % This class implements a transmitter 
    %   Detailed explanation goes here
    
    properties
        no_of_subcarriers;  %Active subcarriers
        total_no_of_subsymbols;
        no_of_subsymbols;
        no_of_blocks;
        no_of_total_subcarriers;
        no_of_bin_data;
        len_of_mapped_signal;
        len_of_CP;
        binary_source;
        maped_signal;
        no_of_bits_per_symbol;
        modulation_type;
        Pulse_Shaping_Filter_Type;
        filter_pulse;
        signal_type;
        gfdm_signal_;
        parallel_signal_;
        transmitted_signal;
        Aii;
        hh;
    end
    
    methods
        function obj = transmitter()
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.no_of_bin_data = 2^16;
            obj.binary_source = randsrc(1,obj.no_of_bin_data,[0 1]);
            obj.no_of_subcarriers = 128;
            obj.no_of_total_subcarriers = 256;
        end
        
        function mapper(obj,mod_type)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here               
            switch mod_type
                case 'QPSK'
                    obj.no_of_bits_per_symbol = 2; 
                    maped_signal_ = qammod(obj.binary_source', 4,'InputType','bit','UnitAveragePower',true);
                    obj.maped_signal = maped_signal_;
                    obj.modulation_type = mod_type;
                    obj.len_of_mapped_signal = length(maped_signal_);
                case '16QAM'
                    obj.no_of_bits_per_symbol = 4;
                    maped_signal_ = qammod(obj.binary_source', 16, 'InputType','bit','UnitAveragePower',true);
                    obj.maped_signal = maped_signal_;
                    obj.modulation_type = mod_type; 
                    obj.len_of_mapped_signal = length(maped_signal_);
            end  
        end
        
        function res = RRCPulseShape(obj)
               M = obj.total_no_of_subsymbols;
               K = obj.no_of_total_subcarriers;
               a = 0.1;
               t = linspace(-M/2, M/2, M*K+1); t = t(1:end-1); t = t';
               g = (sin(pi*t*(1-a))+4*a.*t.*cos(pi*t*(1+a)))./(pi.*t.*(1-(4*a*t).^2));
               g(find(t==0)) = 1-a+4*a/pi;
               g(find(abs(t) == 1/(4*a))) = a/sqrt(2)*((1+2/pi)*sin(pi/(4*a))+(1-2/pi)*cos(pi/(4*a)));

               g = fftshift(g);
               res = g / sqrt(sum(g.*g));
               obj.filter_pulse = res;
        end
        
        function gen_signal = modulate(obj, signal_type_, cp_len, num_of_subsymbols)
            obj.len_of_CP = cp_len;
            maped_sig = obj.maped_signal;
            obj.signal_type = signal_type_;
            if (signal_type_ == 'OFDM')
                obj.no_of_subsymbols = 1;
                obj.no_of_blocks = obj.len_of_mapped_signal/obj.no_of_subcarriers;
                obj.Pulse_Shaping_Filter_Type = '';
                reshaped_signal = reshape(maped_sig, [],obj.no_of_blocks);
                ofdm_signal = [zeros(50,obj.no_of_blocks);reshaped_signal;zeros(78,obj.no_of_blocks)];
                ofdm_signal = ifft(ofdm_signal);
                ofdm_signal_cp = [ofdm_signal(:,end - obj.len_of_CP+1 : end), ofdm_signal];
                obj.hh = ofdm_signal_cp;
                obj.transmitted_signal = reshape(ofdm_signal_cp,1,[]);
                obj.transmitted_signal = obj.transmitted_signal/std(obj.transmitted_signal);
                gen_signal = obj.transmitted_signal;
            end
            if(signal_type_ == 'GFDM')
                obj.Pulse_Shaping_Filter_Type = 'Raised Root Cosine';
                obj.no_of_subsymbols = num_of_subsymbols;
                obj.no_of_blocks = obj.len_of_mapped_signal/(obj.no_of_subcarriers*obj.no_of_subsymbols);
                block_len = obj.no_of_subsymbols * obj.no_of_subcarriers;
                obj.total_no_of_subsymbols = obj.no_of_subsymbols;
                total_block_len = obj.no_of_total_subcarriers * obj.total_no_of_subsymbols;
                GFDM_signal = zeros((obj.no_of_blocks*total_block_len)+cp_len,1);
                for i = 0:obj.no_of_blocks-1
                    sig = maped_sig(block_len*i+1:end-(block_len*(obj.no_of_blocks - 1 - i)));
                    sig_reshaped = reshape(sig, obj.no_of_subcarriers, obj.no_of_subsymbols);
                    sig_total = zeros(obj.no_of_total_subcarriers,obj.total_no_of_subsymbols);
                    sig_total(51:obj.no_of_subcarriers+50,:) = sig_reshaped; 
                    filter_value = RRCPulseShape(obj);
                    sig_rep = repmat(obj.no_of_total_subcarriers*ifft(sig_total),obj.total_no_of_subsymbols,1);
                    sig_sum = zeros(total_block_len, 1);
                    for j = 1:obj.total_no_of_subsymbols
                        sym = sig_rep(:,j).*filter_value;
                        sym = circshift(sym, obj.no_of_total_subcarriers*(j-1));
                        sig_sum = sig_sum + sym;
                    end
                    sig_cp = [sig_sum(end-cp_len+1:end); sig_sum];
               
                    GFDM_signal(i*(total_block_len+cp_len) + (1:total_block_len + cp_len)) = sig_cp;         
                end
                obj.transmitted_signal = GFDM_signal';
                obj.transmitted_signal = obj.transmitted_signal/std(obj.transmitted_signal);
                gen_signal = obj.transmitted_signal;
            end
        end
    end
end

