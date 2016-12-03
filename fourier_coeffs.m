 #Given a function defined on [x0,xf], compute the dim-th partial sum of the Fourrier series expansion of f 
#Method: We are essentially computing the coordinate representation of the vector relative to #the the 2*dim trignometric monomials {cos(pi*x), cos(2*pi*x), cos(3*pi*x), ..., cos(n*pi*x), sin(pi*x), sin(2*pi*x), sin(3*pi*x),..., sin(n*pi*x) }, which span a function space of dimension n (L2).

function [ccoord, scoord]  = fourier_coeffs( dim, x0, xf, fx ) 


	#the interval of integration is [x0,xf]. The fineness of the partition (h=coarse) can be changed here, if needed
  coarse=0.001;
  xx=[x0:coarse:xf];
  nn=length(xx);
  
  #ccoord is the container of coefficients the to each cosine wave
  #scoord is the container of coefficients the to each sine wave
  ccoord=scoord=zeros(1,dim+1);   
  
  ccoord(1)= trapz(xx, fx(xx)  ) / pi;
  scoord(1)=0; 
  
  for i=2:dim+1
	  #cos_n contains the values of each cosine wave
    #sin_n contains the values of each sine wave
	  cos_n= cos( (i-1)*xx );
	  sin_n=  sin( (i-1)*xx );
	  
    #Compute the integrands before integration  	
	  temp_f1 = fx(xx) .* cos_n;
		temp_f2 = fx(xx) .* sin_n;
	    	
    #Compute the i-th coefficients by integration   	    	  
		ccoord(i)=  trapz(xx,temp_f1) / pi;
		scoord(i)=  trapz(xx,temp_f2) / pi;

	end



end
