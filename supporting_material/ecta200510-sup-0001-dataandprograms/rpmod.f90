module globals
	use myfuncs
	use dotops
	
	integer, parameter	:: ngq=40
	real	:: GQxw(ngq,2)
	
	real				:: level=0.05
	integer, parameter	:: qmax=30, nc=30

	end module
	
module rpmod
	use globals
	implicit none
	
	real, allocatable	:: global_Omcv(:,:)
	real				:: global_level
	
	contains
	
	function getrp_fromev(ev) result(val)
		real	:: ev(:),val
		real	:: u
		integer	:: i,q
		q=size(ev)
		if(q==0) then
			val=-100
			return
		endif
		val=0
		do i=1,ngq
			u=GQxw(i,1)
			val=val+GQxw(i,2)/sqrt((1-u**2)*product(1.0+ev/(1-u**2)))
		enddo
		val=2*val/pi
	end function
	
	function getevs(Om,cv) result(val)
		use evlrg_int
		real	:: Om(:,:),cv
		real, allocatable	:: val(:)
		real	:: A(size(Om,1),size(Om,1))
		complex	:: evsc(size(Om,1))
		A(:,1)=Om(:,1)
		A(:,2:)=-cv**2*Om(:,2:)
		call evlrg(A,evsc)
		val=-pack(real(evsc),real(evsc)<0)/maxval(real(evsc))	
	end function
	
	function getrp(Om,cv) result(val)
		real	:: Om(:,:),cv,val
		val=getrp_fromev(getevs(Om,cv))
	end function
	
	function getStudentcv(q,level) result(val)
		integer	:: q
		real	:: level,val
		val=sqrt(tin(1-level/2,real(q))**2/q)
	end function

	function getcv(Om,level) result(val)
		use zuni_int
		real	:: Om(:,:),level, val
		real	:: temp
		real, external	:: ext_getlcv
		
		global_Omcv=Om
		global_level=level
		temp=log(getStudentcv(size(Om,1)-1,level))-log(2.0)
		val=temp+4*log(2.0)
		call zuni(ext_getlcv,temp,val,TOL=0.001)	
		val=exp(val)
	end function
	
	function getlength_fromev(ev) result(val)
		real	:: ev(:),val
		integer	:: j
		real	:: s,t
		val=0
		do j=1,ngq
			s=GQxw(j,1)
			t=(1-1/s)**2
			val=val+GQxw(j,2)*((1.0-s)/s**3)*sum(ev/(1+2*t*ev))/sqrt(t*product(1+2*t*ev))
		enddo
		val=4*val/sqrt(pi)
	end function
		
	function getlength(Om,cv) result(val)
		use evlsf_int
		real	:: Om(:,:),cv,val
		real	:: ev(size(Om,1)-1)
		
		call evlsf(Om(2:,2:),ev)
		ev=ev*cv**2/Om(1,1)
		val=getlength_fromev(ev)
	end function

	function getIlength(q,cv) result(val)
		integer	:: q
		real	:: cv,val
		real	:: ev(q)
		ev=cv**2
		val=getlength_fromev(ev)
	end function
		
end module
	
function ext_getlcv(lcv) result(val)
	use rpmod
	implicit none
	real	:: lcv,val
	val=getrp(global_Omcv,exp(lcv))-global_level
end function
	
	