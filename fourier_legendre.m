#Given a a function defined on a real interval [x0,xf],  determine the coordinate representation of the function relative to the n-dimensional #Legendre basis
#Method: we are essentially computing the dim-th partial sum of the Fourrier-Legendre series expansion of f
#	 That expansion is done relative the n first Legendre polynomials, which are ONLY defined on [-1,1]. If anyone wants to help with scaling,   please fork.	
function coord = fourier_legendre( x0, xf, dim, f )
  
	#the interval of integration is [x0,xf]. The fineness of the partition (h) can be changed here, if needed
  	coarse=0.001;
  	xx=[x0:coarse:xf];
  	nn=length(xx);
  
  	#coord is the container of the vector of coefficients to each legendre polynomial in the expansion 
  	coord=zeros(1,dim+1);   
  
	for i=1:dim+1
	    	#We first compute the ith Legendre polynomial of order 0
	    	temp= (legendre( i - 1,xx))(1,:)  ;  
	    	#Compute the values of  f(x)*Pn(x) and Pn(x)^2
	    	tempf=f(xx).*temp;
	    	temp2=temp.*temp;
	    
	    	#Compute the i-th coefficient by integration   
	    	#formula: integral( f(x)*Pn(x) ) / integral( Pn(x)^2 )   
		coord(i)=  trapz(xx,tempf) / trapz(xx,temp2);
	end

  
endfunction




