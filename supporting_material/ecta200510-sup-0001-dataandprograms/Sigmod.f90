module Sigmod
	use globals
	implicit none
	
	contains
	
	function getOmfromW(W,Sig) result(val)
		real	:: W(:,:),Sig(:,:),val(size(W,2),size(W,2))
		val=matmul(transpose(W),matmul(Sig,W))
	end function
	
	function getcSig(c,s) result(val)
		real	:: c,s(:,:)
		real, allocatable	:: val(:,:)
		integer	:: n,i,j
		n=size(s,2)
		allocate(val(n,n))
!$omp parallel do private(j)
		do i=1,n
			do j=1,n
				val(j,i)=exp(-c*norm2(s(:,i)-s(:,j)))
			enddo
		enddo
	end function	

	function getcpSig(c,p,s) result(val)
		real	:: c,s(:,:)
		integer	:: p
		real, allocatable	:: val(:,:)
		real	:: d
		integer	:: n,i,j
		n=size(s,2)
		allocate(val(n,n))
!$omp parallel do private(j,d)
		do i=1,n
			do j=1,n
				d=norm2(s(:,i)-s(:,j))
				select case(p)
					case(0)
						val(j,i)=exp(-c*d)
					case(1)
						val(j,i)=(1+c*d)*exp(-c*d)
					case(2)
						val(j,i)=(1+c*d+(c*d)**2/3.0)*exp(-c*d)
					case default
						val(j,i)=exp(-c*d**2)
				end select
			enddo
		enddo
	end function	

	function getsmax(s) result(smax)
		real	:: s(:,:),smax
		integer	:: l1,l2,n
		n=size(s,2)
		smax=0
		do l1=1,n
			do l2=l1+1,n
				smax=max(smax,norm2(s(:,l1)-s(:,l2)))
			enddo
		enddo
	end function

	function getWeigen(Sig,q) result(val)
		use evesf_int
		real	:: Sig(:,:)
		integer	:: q
		real	:: c,val(size(Sig,1),q+1)
		real	:: Sigm(size(Sig,1),size(Sig,1)),evs(q)
		integer	:: i,n
		n=size(Sig,1)
		do i=1,n
			Sigm(:,i)=Sig(:,i)-sum(Sig(:,i))/n
		enddo
		do i=1,n
			Sigm(i,:)=Sigm(i,:)-sum(Sigm(i,:))/n
		enddo
		call evesf(q,Sigm,.false.,evs,val(:,2:q+1))
		val(:,1)=1.0/sqrt(real(n))
	end function
	
	function getKVSig(s,b) result(Sig)
		real	:: s(:,:),b,Sig(size(s,2),size(s,2))
		integer	:: i,j,n
		n=size(s,2)
		do i=1,n
			do j=1,n
				Sig(j,i)=max(1.0-norm2(s(:,i)-s(:,j))/b,0.0)
			enddo
		enddo
		do i=1,n
			Sig(:,i)=Sig(:,i)-sum(Sig(:,i))/n
		enddo
		do i=1,n
			Sig(i,:)=Sig(i,:)-sum(Sig(i,:))/n
		enddo
	end function
	
	function getWKV(s,b) result(val)
		use evesf_int
		integer, parameter	:: qmax=200
		real	:: s(:,:),b
		real, allocatable	:: val(:,:), Sig(:,:)
		real	:: evs(qmax), evcs(size(s,2),qmax),rt
		integer	:: i,j,n,q
		n=size(s,2)
		Sig=getKVSig(s,b)
		call evesf(qmax,Sig,.false.,evs,evcs)
		evs=evs/sum(evs,evs>0)
		rt=0;q=0
		do i=1,qmax
			rt=rt+evs(i)
			q=q+1
			if(rt>0.995) exit
		enddo
		allocate(val(n,q+1))
		evs=evs/(sum(evs(1:q))/q)
		do i=1,q
			val(:,i+1)=.2*sqrt(evs(i))*evcs(:,i)
		enddo
		val(:,1)=1.0/sqrt(real(n))
!		print *,"smallest eigenvalue of KV spectral decomposition", evs(q)/evs(1)
	end function
	
	function getWSun(s,k1,k2) result(val)
		real	:: s(:,:)
		integer	:: k1,k2
		real, allocatable	:: val(:,:)
		integer	:: i1,i2,q,n,k1c,k2c
		real	:: st(size(s,1),size(s,2))
		n=size(s,2)
		st(1,:)=s(1,:)/maxval(s(1,:))
		st(2,:)=s(2,:)/maxval(s(2,:))
		
		allocate(val(n,1))
		val=1
		do i1=0,k1
			do i2=0,k2
				if(i1==0 .and. i2==0) cycle
				val=reshape([val,cos(2*Pi*matmul(real([i1,i2]),st)),sin(2*Pi*matmul(real([i1,i2]),st))],[n,size(val,2)+2])
			enddo
		enddo
		val=matmul(val,invertgen(transpose(choleski(matmul(transpose(val),val)))))
	end function
	
	function getavcorr(Sig) result(val)
		real	:: Sig(:,:),val
		integer	:: i,j,n
		n=size(sig,1)
		val=0
		do i=1,n
			do j=1,i-1
				val=val+Sig(j,i)
			enddo
		enddo
		val=val/(n*(n-1)/2)
	end function
	
	function getavc_c(avc,p,s) result(c)
		use rootfinder
		real	:: avc,s(:,:),c
		integer	:: p
		integer	:: n,i,j,l
		real	:: c0,c1,tol,f
		real, allocatable	:: dns(:)
		integer(4)	:: rstat
		
		n=size(s,2)
		allocate(dns(n*(n-1)/2))
		l=1
		do i=1,n
			do j=1,i-1
				dns(l)=norm2(s(:,i)-s(:,j))
				l=l+1
			enddo
		enddo
		c0=10; c1=80
		do
			if(getavcd(c0)*getavcd(c1)<0) exit
			c0=c0/2; c1=c1*2
		enddo
		rstat=0; tol=0.001
		call zero_rc(c0,c1,tol,c,rstat,f)
		do
			f=getavcd(c)
			call zero_rc(c0,c1,tol,c,rstat,f)
			if(rstat==0) exit
		enddo
		
	contains
		function getavcd(c) result(val)
			real	:: c,val
			select case(p)
				case(0)
					val=sum(exp(-c*dns))/size(dns)
				case(1)
					val=sum((1+c*dns)*exp(-c*dns))/size(dns)
				case(2)
					val=sum((1+c*dns+(c*dns)**2/3.0)*exp(-c*dns))/size(dns)
				case default
					val=sum(exp(-c*dns**2))/size(dns)
			end select
			val=val-avc	
		end function
	end function	

end module
			
