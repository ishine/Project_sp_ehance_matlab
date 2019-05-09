fin = '../voice/omtest2'; %040_ANC_PC2.wav

[y,Fs]=audioread([fin '.wav']);  % read size of input data, Fs and NBITS

out = omlsa_func (y,Fs); 


 audiowrite([[fin,'omlsa_cal_out.wav'] '.wav'],out,16000);