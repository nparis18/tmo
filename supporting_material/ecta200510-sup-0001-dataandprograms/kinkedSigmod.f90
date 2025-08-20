module kinkedSigmod
	use Sigmod
	use compute
	implicit none
	
	real	:: global_shift, global_c02, global_r, global_c12
!$omp threadprivate(global_shift, global_c02, global_r, global_c12)
	
	contains
	
	function getSdelcorr(r,lam0,c0,initflag) result(val)
		use qdag_int
		use qdval_int
		real	:: r,lam0,c0,val
		logical				:: initflag
		integer				:: i
		integer, parameter	:: nr=100
		real, parameter		:: rgrid(nr)=[(sqrt(2.0)*((i-1.0)/(nr-1))**2,i=1,nr)]
		real, save			:: fgrid(nr)
		real, external		:: ext_Sdelcorr
		if(initflag) then
			global_shift=specdens(lam0,c0**2)
			global_c02=c0**2
			do i=1,nr
				global_r=rgrid(i)
				call qdag(ext_Sdelcorr,0.0,lam0,fgrid(i),errabs=0.0001,irule=6)
			enddo
			val=-1
			return
		endif
		val=qdval(r,rgrid,fgrid,check=.false.)
	end function
	
	function getSdelSig(lam0,c0,s,Sig0) result(val)
		use evasf_int
		real	:: lam0,c0,s(:,:),Sig0(:,:)
		real, allocatable	:: val(:,:)
		real	:: evdiff(1)
		integer	:: n,i,j
		n=size(s,2)
		allocate(val(n,n))
		val(1,1)=getSdelcorr(0.0,lam0,c0,.true.)
!$omp parallel do private(j) schedule(dynamic)
		do i=1,n
			do j=1,i
				val(j,i)=getSdelcorr(norm2(s(:,i)-s(:,j)),lam0,c0,.false.)
				val(i,j)=val(j,i)
			enddo
		enddo
		call evasf(1,val-Sig0/Sig0(1,1),small=.true.,eval=evdiff)
		do i=1,n
			val(i,i)=val(i,i)-evdiff(1)
		enddo
	end function
	
	function getc1delcorr(r,lam0,c0,c1,initflag) result(val) ! generates kinked spectral densities
		use qdag_int
		use qdval_int
		real	:: r,lam0,c0,c1,val
		logical				:: initflag
		integer				:: i
		integer, parameter	:: nr=100
		real, parameter		:: rgrid(nr)=[(sqrt(2.0)*((i-1.0)/(nr-1))**2,i=1,nr)]
		real, save			:: fgrid(nr)
		real, external		:: ext_c1delcorr
		if(initflag) then
			global_shift=specdens(lam0,c0**2)/specdens(lam0,c1**2)
			global_c02=c0**2
			global_c12=c1**2
			do i=1,nr
				global_r=rgrid(i)
				call qdag(ext_c1delcorr,0.0,lam0,fgrid(i),errabs=0.0001,irule=6)
			enddo
			val=-1
			return
		endif
		val=qdval(r,rgrid,fgrid,check=.false.)
	end function
	
	function getc1Sig(lam0,c0,c1,s) result(val)
		real	:: lam0,c0,c1,s(:,:)
		real, allocatable	:: val(:,:),Sig1(:,:)
		integer	:: n,i,j
		n=size(s,2)
		allocate(val(n,n))
		val(1,1)=getc1delcorr(0.0,lam0,c0,c1,.true.)
!$omp parallel do private(j) schedule(dynamic)
		do i=1,n
			do j=1,i
				val(j,i)=getc1delcorr(norm2(s(:,i)-s(:,j)),lam0,c0,c1,.false.)
				val(i,j)=val(j,i)
			enddo
		enddo
		Sig1=getcSig(c1,s)
		val=val+global_shift*Sig1/Sig1(1,1)
	end function


	function specdens(lam,c2) result(val)
		real	:: lam,c2,val
		val=1.0/(c2*(1.0 + lam**2/c2)**(1.5))
	end function
	
end module

function ext_sdelcorr(lam) result(val)
	use kinkedSigmod
	real	:: lam,val
	val=lam*Bessel_J0(lam*global_r)*(specdens(lam,global_c02)-global_shift)
end function

function ext_c1delcorr(lam) result(val)
	use kinkedSigmod
	real	:: lam,val
	val=lam*Bessel_J0(lam*global_r)*(specdens(lam,global_c02)-global_shift*specdens(lam,global_c12))
end function
	
	