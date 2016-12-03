#Author: Koffi Dossa
#Given a function defined on [x0,xf], return a function handle of the dim-th partial sum of the Fourrier series expansion of f 
#The closed form solution of the fourier series coefficients are well known. See the litterature
function fourier_handle  = fourier_exp( dim, x0, xf, fx ) 


  #the interval of integration is [x0,xf]. The fineness of the partition (h=coarse) can be changed here, if needed
  coarse=0.001;
  xx=[x0:coarse:xf];
  nn=length(xx);
  
  #ccoord is the container of coefficients the to each cosine wave
  #scoord is the container of coefficients the to each sine wave
  ccoord=scoord=zeros(1,dim+1);   
  
  #First compute the constant term
  ccoord(1)= trapz(xx, fx(xx)  ) / pi;
  scoord(1)=0; 
  
  #You should know what an inner product is in order to fully appreciate what is going on in this loop
  for i=2:dim+1
	#cos_n contains the values of each cosine wave
    	#sin_n contains the values of each sine wave
	cos_n= cos( (i-1)*xx );
	sin_n=  sin( (i-1)*xx );
	  
    	#Compute the integrands before integration  	
	temp_f1 = fx(xx) .* cos_n;
	temp_f2 = fx(xx) .* sin_n;
	    	
    	#Compute the i-th coefficients by integration   	    	  
	ccoord(i)=  trapz(xx,temp_f1) / pi;;
	scoord(i)=  trapz(xx,temp_f2) / pi;

	end

  fourier_handle = @(x) ccoord(1);
  
  #Constructing the function (as a sum)
  for i=2:dim+1    
    fourier_handle= @(x) fourier_handle(x) + ccoord(i)*cos( (i-1)*x) + scoord(i)*sin((i-1)*x);
  end
  
end
