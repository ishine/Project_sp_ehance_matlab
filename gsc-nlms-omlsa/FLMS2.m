function   FLMS2(  ) 
 
% papare the input data
finl = './voice/t15_U_B';
fout = [finl '_flms2_div12'];
 
[x,fs]= audioread([finl '.wav']);  % main mic
%[x2,fs2]= audioread([finr '.wav']); % ref mic
lambda = 0.88;
alpha = 0.98/12;
filter_length = 512;
input_signal = x(:,2);
desired_signal = x(:,1);


% initialization 
FILTER_COEFF = zeros(2*filter_length,1); 
input_length = length(input_signal); 
estimated_power = zeros(2*filter_length,1); 
 
%length of signal should be integer multiple of filter_length 
block_length = floor(input_length/filter_length)*filter_length; 
 
%truncate signal 
input_signal = input_signal(1:block_length); 
desired_signal = desired_signal(1:block_length); 
 
% make sure that input signal and desired signal are column vectors 
input_signal = input_signal(:); 
desired_signal = desired_signal(:); 
 
error = desired_signal; 
 
% no.of blocks 
Blocks = block_length/filter_length; 
 
out_y = zeros(input_length,1);

% omlsa init
N_OMLSA = 512;
q_OMLSA = N_OMLSA/4;
Y_om_in = zeros(N_OMLSA,1);
U_om_in = zeros(N_OMLSA,1);
yout_omlsa = zeros(input_length,1);
k_omlsa = 0;

y_omlsa_flms_vec = zeros(2*filter_length,1); 
y_omlsa_flms = zeros(2*filter_length,1);

% power ratio 
pow_om_array = zeros(Blocks,1);
pow_err_array = zeros(Blocks,1);
pow_ratio_array = zeros(Blocks,1);
pow_ratio_num_array= zeros(Blocks,1);
% pow_om= zeros(2*filter_length,1);  pow_err= zeros(2*filter_length,1);
% pow_ratio_num = ones(2*filter_length,1);
pow_om  = 0;
pow_err =0;
pow_ratio_num = 0;
 pow_ratio_num = 1;
beta = 1;


% loop, FLMS 
for k=1:Blocks-1 
 
    % block k-1, k; transformed input signal U(k) 
    INPUT_SIGNAL = fft( input_signal((k-1)*filter_length+1:(k+1)*filter_length) ,2*filter_length); 
 
    % block k, output signal y(k), last M elements 
    output = ifft(INPUT_SIGNAL.*FILTER_COEFF); 
    output = output(filter_length+1:2*filter_length,1); 
 
    % block k; desired signal 
    desired_vec = desired_signal(k*filter_length+1:(k+1)*filter_length); 
     
    % block k, error signal 
    error(k*filter_length+1:(k+1)*filter_length,1) = desired_vec-output; 
    
    out_y(k*filter_length+1:(k+1)*filter_length,1) = output;
 
    % transformation of estimation error 
    ERROR_VEC = fft([zeros(filter_length,1);error(k*filter_length+1:(k+1)*filter_length)],2*filter_length); 
      
    p_mode = 2;
     
    if(p_mode ==1)
    % estimated power 
    % 1£© orig pest
     estimated_power=lambda*estimated_power+(1-lambda)*abs(INPUT_SIGNAL).^2; 
    else
    % 2) combined omlsa pest  
       estimated_power=lambda*estimated_power+(1-lambda)*abs(INPUT_SIGNAL).^2; 
     
         beta = pow_ratio_num;
       if( pow_ratio_num < 0.4)
         beta  = 0;% pow_ratio_num*0.5;
       end
   
    end
    % block k, inverse of power 
    DESIRED_VEC = 1./(1e-2+estimated_power) ; 
    
   %  DESIRED_VEC = estimated_power.^2 ./ (1e-2+ estimated_power.^2);
 
    % estimated gradient 
    phivec = ifft(DESIRED_VEC.*conj(INPUT_SIGNAL).*ERROR_VEC,2*filter_length); 
    phivec = phivec(1:filter_length); 
 
    % update of weights 
    FILTER_COEFF = FILTER_COEFF+ beta * alpha * fft([phivec;zeros(filter_length,1)],2*filter_length); 
 
    % The error vector should have only real values. 
    error = real(error(:)); 
    
    % prepare the input for omlsa
    err_temp = error(k*filter_length+1:(k+1)*filter_length);
    x_new    = input_signal(k*filter_length+1:(k+1)*filter_length);
    
        for tmp = 1:4   % 4 *128                            
            Y_om_in  = [ Y_om_in(q_OMLSA+1:N_OMLSA); err_temp((tmp-1)*q_OMLSA+1:tmp*q_OMLSA)];     
            U_om_in  = [ U_om_in(q_OMLSA+1:N_OMLSA); x_new((tmp-1)*q_OMLSA+1:tmp*q_OMLSA)];
 
            yout_omlsa_t  =  flms_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );      
            yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;   k_omlsa = k_omlsa+1;            
            
             
           y_omlsa_flms = [  y_omlsa_flms(q_OMLSA+1 : length(y_omlsa_flms)  ) ;yout_omlsa_t];        

        end   
        
           
        % cal the power of omlsa output,  ratio between  error omlsa 
        y_omlsa_flms_vec = fft(y_omlsa_flms , 2*filter_length);
        
        
%            pow_om =  0.88* pow_om  +   0.12 *abs(y_omlsa_flms_vec).^2;               
%            pow_err = 0.88* pow_err   +   0.12 * abs(ERROR_VEC).^2;    
           pow_om   =  0.88 * pow_om  +   0.12 * (y_omlsa_flms_vec' * y_omlsa_flms_vec) ;               
          pow_err =  0.88 * pow_err  +   0.12 * abs(ERROR_VEC' * ERROR_VEC) ;      

        pow_ratio = pow_err ./ pow_om;
         
        pow_ratio_num = pow_ratio ./ (pow_ratio + 2 );
         
        if(k>3)
        pow_om_array(k)    = pow_om;
        pow_err_array(k)   = pow_err;
        
        
        pow_ratio_array(k) = pow_ratio;
       
        
        pow_ratio_num_array(k) = pow_ratio_num;        
         end
end
 
D_E = [  desired_signal,error];
audiowrite([fout '_D_E.wav'],D_E,fs);

D_Y = [  desired_signal,out_y(1:length(desired_signal))]; 
audiowrite([fout '_D_Y.wav'],D_Y,fs);

E_B = [ error, input_signal];
audiowrite([fout '_E_B.wav'],E_B,fs);

audiowrite([fout '_om.wav'],yout_omlsa,fs);






% figure
%    figure(1)
%    x = 1:Blocks  ;
%    subplot(411);  plot(x, pow_err_array, 'r'); title('1) error power')
%    subplot(412);  plot(x, pow_om_array, 'r'); title('2) omlsa power')
%    subplot(413);  plot(x, pow_ratio_array, 'r'); title('3) error om  power ratio')
%    subplot(414);  plot(x, pow_ratio_num_array, 'r'); title('4) power ratio num')


fprintf('end of flms2\n');
