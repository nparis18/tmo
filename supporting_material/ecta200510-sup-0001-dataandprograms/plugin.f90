module plugin
	use optmod
	implicit none

	integer, parameter	:: np=3,nc=30
	
	real, allocatable	:: iSigs(:,:,:,:), iSigvecs(:,:,:)
	real				:: Sigprops(2,nc,np)
	
	private np,nc,iSigs,Sigprops,iSigvecs

	contains
	
	

	subroutine setSigs(s,avc0)
		real	:: s(:,:),avc0
		real, allocatable	:: Sig(:,:)
		real	:: avc,c0,c
		integer	:: n,ip,ic
		
		n=size(s,2)
		if(allocated(iSigs)) deallocate(iSigs,iSigvecs)
		allocate(iSigs(n,n,nc,np),iSigvecs(n,nc,np))
		do ip=1,np
			c0=getavc_c(avc0,ip-1,s)
!$omp parallel do private (c,Sig)
			do ic=1,nc
				c=c0*exp(-1+4.0*ic/real(nc))
				Sig=getcpSig(c,ip-1,s)				
				Sigprops(1,ic,ip)=logdetpd(Sig)
				Sigprops(2,ic,ip)=sum(Sig)
				call linds(Sig,iSigs(:,:,ic,ip))
				if(iercd()>0) then
					iSigs(:,:,ic,ip)=eye(n)
					Sigprops(1,ic,ip)=1E50
					iSigvecs(:,ic,ip)=0
					cycle
				endif
				iSigvecs(:,ic,ip)=sum(iSigs(:,:,ic,ip),dim=1)/sum(iSigs(:,:,ic,ip))
			enddo
		enddo
	end subroutine
	
	function getMLErp(s,avc) result(val)
		real	:: s(:,:),avc,val
		integer, parameter	:: ndgp=5,nsim=20000
		real	:: c0,lam0,y(size(s,2)),rp(2,ndgp),Sig0(size(s,2),size(s,2)),chol(size(s,2),size(s,2))
		real	:: s2,sy
		integer	:: i,l,n
		
		n=size(s,2)
		c0=getavc_c(avc,0,s)
		lam0=5*c0
		Sig0=getcSig(c0,s)
		call setSigs(s,avc)
		rp=0
		do i=1,ndgp
			chol=getc1Sig((real(i)/ndgp)*lam0,c0,4*c0,s)
			chol=chol/(sum(chol)/real(n))
			chol=choleski(chol)
!$omp parallel do private(y,s2,sy) reduction (+:rp)			
			do l=1,nsim
				call rnnoa(y)
				y=matmul(chol,y)
				sy=sum(y)
				s2=s2MLE(y)
				rp(2,i)=rp(2,i)+sqrt(s2/n)
				if(sy**2>1.96**2*s2) rp(1,i)=rp(1,i)+1
			enddo
		enddo
		rp=rp/nsim
		val=maxval(rp(1,:))	
	end function
	
	function s2MLE(y) result(val)
		real	:: y(:),val
		integer	:: ip,ic,i(2)
		real	:: ll(nc,np)
		do ip=1,np
			do ic=1,nc
				y=y-sum(iSigvecs(:,ic,ip)*y)
				ll(ic,ip)=-.5*size(y)*log(sum(y*matmul(iSigs(:,:,ic,ip),y)))-.5*Sigprops(1,ic,ip)
			enddo
		enddo
		i=maxloc(ll)
		if(i(1)==1 .or. i(1)==nc) then
			print *,i
		endif
		y=y-sum(iSigvecs(:,i(1),i(2))*y)
		val=(sum(y*matmul(iSigs(:,:,i(1),i(2)),y))/size(y))*Sigprops(2,i(1),i(2))
	end function	
	
end module


	