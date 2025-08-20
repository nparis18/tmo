module optmod
	use kinkedSigmod
	use compute
	implicit none
	
	integer, parameter	:: qx=100, n0=20
	real		:: lam0grid(n0)
	
	contains
	
	function s2opt(s,c0,c1,lam0) result(val)		! c1<0 indicates upper bound on expected length under white noise 
		real	:: s(:,:), c0,c1,lam0,val(2)
		integer, parameter	:: n0=20,nsim=50000,nds=20
		real :: dgrid(nds)
		real, allocatable	:: Sig0(:,:),Sig1(:,:),W(:,:),Om(:,:)
		integer	:: q
		real	:: cv
		real	:: condm(qx,n0),conds2(n0),chols(qx,qx,n0),Omis(qx,qx,0:n0),Omdets(0:n0),Y(qx)
		real, allocatable	:: Ys(:,:),LR(:,:),lagr(:,:,:)
		real	:: fprop,la(nds),mux,sx,p
		real	:: rpt(n0),lam(n0),rp(n0),d,avl
		integer		:: i,l,j,i2,ic

		Sig0=getcSig(c0,s)
		W=getWeigen(Sig0,30)
		call setSCPCqcv(W,s,c0,q,cv)	

		W=getWeigen(Sig0,qx)

		do i=1,n0
			lam0grid(i)=(real(i)/n0)*lam0
			if(c1>0) then
				Om=getOmfromW(W,getc1Sig(lam0grid(i),c0,c1,s))
			else
				Om=getOmfromW(W,getSdelSig(lam0grid(i),c0,s,Sig0))
			endif
			Om=Om/Om(1,1)
			Omis(:,:,i)=invertpd(Om(2:,2:))
			condm(:,i)=matmul(Omis(:,:,i),Om(2:,1))
			conds2(i)=Om(1,1)-sum(condm(:,i)*Om(2:,1))
			chols(:,:,i)=choleski(Om(2:,2:))
			Omdets(i)=detpd(Om(2:,2:))
			rpt(i)=max(getrp(Om(1:q+1,1:q+1),cv),level)
		enddo
		if(c1>0) then
			Om=getOmfromW(W,getcSig(c1,s))
			Om=Om/Om(1,1)
		else
			Om=eye(qx+1)
		endif
		Omis(:,:,0)=invertpd(Om(2:,2:))
		Omdets(0)=detpd(Om(2:,2:))
		
		allocate(LR(1:n0,nsim),lagr(nds,0:n0,nsim))
		do i=1,nds
			dgrid(i)=(4.5*(i-1))/nds
		enddo
		
!$omp parallel do private(Y,i,mux,sx,j,fprop)		
		do l=1,nsim
			call rnnoa(Y)
			Y=matmul(chols(:,:,mod(l,n0)+1),Y)
			Y=sqrt(real(q))*Y/norm2(Y)
			sx=sum(Y*matmul(Omis(:,:,0),Y))
			lagr(:,0,l)=dgrid*sqrt(real(qx))*(sx)**(-.5*(qx+1))/sqrt(Omdets(0))
			do i=1,n0
				sx=sum(Y*matmul(Omis(:,:,i),Y))
				LR(i,l)=(sx)**(-.5*qx)/sqrt(Omdets(i))
				mux=sum(Y*condm(:,i))
				sx=sqrt(sx*conds2(i)/qx)
				do j=1,nds
					lagr(j,i,l)=1-gausscdf((dgrid(j)-mux)/sx)+gausscdf((-dgrid(j)-mux)/sx)
				enddo
			enddo
			fprop=sum(LR(1:n0,l))/n0
			LR(:,l)=LR(:,l)/fprop
			lagr(:,0,l)=lagr(:,0,l)/fprop
			do i=1,n0
				lagr(:,i,l)=lagr(:,i,l)*LR(i,l)
			enddo
		enddo
		lam=1
		rp=level
		do ic=1,2000
			lam=lam*exp(2*(rp-rpt))
			rp=0
!$omp parallel do private(la,i,d,i2,p) reduction(+:rp)		
			do l=1,nsim
				la=lagr(:,0,l)+matmul(lagr(:,1:,l),lam)
				i=minloc(la,dim=1)
				if(i>1 .and. i<nds) then
					d=getquadmin(dgrid(i-1:i+1),la(i-1:i+1))
					if(d>dgrid(i)) then
						i2=i+1
					else 
						i2=i-1
					endif
					p=abs((d-dgrid(i2))/(dgrid(i2)-dgrid(i)))
				else
					if(ic==2000 .and. i==nds) print *,"x"
					i2=i
					p=1
				endif				
				rp=rp+p*lagr(i,1:,l)+(1-p)*lagr(i2,1:,l)
			enddo
			rp=rp/nsim

			if(mod(ic,2000)==-1) then
				call mdisp(lam0grid)
				call mdisp(rpt.cvr.rp.cud.lam)
			endif
		enddo
		avl=0
		do l=1,nsim
			la=lagr(:,0,l)+matmul(lagr(:,1:,l),lam)
			i=minloc(la,dim=1)
			if(i>1 .and. i<nds) then
				d=getquadmin(dgrid(i-1:i+1),la(i-1:i+1))
				if(d>dgrid(i)) then
					i2=i+1
				else 
					i2=i-1
				endif
				p=abs((d-dgrid(i2))/(dgrid(i2)-dgrid(i)))
			else
				i2=i
				p=1
			endif
			avl=avl+p*lagr(i,0,l)+(1-p)*lagr(i2,0,l)
		enddo
		avl=2*avl/nsim
		val=[avl,maxval(rpt)]
	end function

	
end module

	