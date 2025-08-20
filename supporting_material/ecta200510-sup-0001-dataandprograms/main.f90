module runmod
!dir$ nooptimize
	use compute
	use robmod
	use optmod
	implicit none
	
	integer						:: n
	character(len=100), parameter	:: dir="C:/Dropbox/mystuff/SpatialHAR/play/"
	integer, parameter			:: ndraws=5

	contains
	
	function getclusterW(s,nc) result(val)
		real	::s(2,n),val(n,nc+1)
		integer	:: nc
		real	:: s_c(2,nc), sdist(nc,n), cdist(nc)
		
		real	:: smin(2),smax(2)
		integer	:: nccounts(nc),cc,mc,j,l,c,lc
		
		if(.not.(nc==4 .or. nc==9)) then
			print *,"value of nc in getcluster not supported"
			stop
		endif
		smin=minval(s,dim=2); smax=maxval(s,dim=2)
		s_c(:,1)=smin
		s_c(:,2)=smax
		s_c(:,3)=[smin(1),smax(2)]
		s_c(:,4)=[smax(1),smin(2)]
		if(nc==9)then
			s_c(:,5)=[.5*smin(1)+.5*smax(1),smin(2)]
			s_c(:,6)=[.5*smin(1)+.5*smax(1),smax(2)]
			s_c(:,7)=[smin(1),.5*smin(2)+.5*smax(2)]
			s_c(:,8)=[smax(1),.5*smin(2)+.5*smax(2)]
			s_c(:,9)=[.5*smin(1)+.5*smax(1),.5*smin(2)+.5*smax(2)]
		endif
		do l=1,n
			do j=1,nc
				sdist(j,l)=norm2(s(:,l)-s_c(:,j))
			enddo
		enddo
		val=0; nccounts=0; cc=0; mc=mod(n,nc)
		do l=1,n
			do j=1,nc
				cdist(j)=minval(sdist(j,:),dim=1)
			enddo
			j=maxloc(cdist,dim=1)
			lc=minloc(sdist(j,:),dim=1)
			val(lc,j+1)=1
			sdist(:,lc)=1E10
			nccounts(j)=nccounts(j)+1
			if(nccounts(j)>=merge(n/nc+1,n/nc,cc<mc)) then
				sdist(j,:)=-1E10
				cc=cc+1
			endif
		enddo
		do j=2,nc+1
			val(:,j)=val(:,j)/nccounts(j-1)
		enddo
		val(:,1)=sum(val(:,2:),dim=2)/nc
		do j=2,nc+1
			val(:,j)=val(:,j)-val(:,1)
		enddo
		val(:,1)=val(:,1)/(sum(val(:,1)/sqrt(real(n))))
!		call mdisp(val)
	end function

	subroutine locerror
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		logical		:: iflag
		real		:: avrp
		real		:: rp(3)
		real		:: out(size(rp),nstates*ndraws)
		real		:: tab(3,size(rp),2)
		
		avc0=0.03
		do ilp=1,2
			lp=ilp-1
			print *,"ooooooooooooooo locerror oooooooooooooooooooooooooooo"
			print *,"avc0=",avc0," lp=",lp
			out=-1; is=1
			do i=1,nstates
				do id=1,ndraws
					s=getstatelocs(i,lp,n)
					call getmrp
					
					out(:,is)=[rp]
					is=is+1
				enddo
			enddo
			call mdisp(out)
			do i=1,size(rp)
				tab(:,i,ilp)=quantile_v(out(i,:),[.05,.5,.95])
			enddo
			call mdisp(tab(:,:,ilp))
			call printtime
		enddo
		
	contains
		
		subroutine getmrp
			real, parameter	:: vlist(3)=[0.75,1.5,2.25]
			real, allocatable:: Om(:,:), W(:,:), Sig(:,:)
			s0=s
			c00=getavc_c(avc0,0,s0)
			do l=1,3
				call rnun(s(1,:))
				call rnun(s(2,:))
				s=s0+vlist(l)*2*(s-.5)/40
				c0=getavc_c(avc0,0,s)			
				Sig=getcSig(c0,s)
				W=getWeigen(Sig,qmax+1)
				call setSCPCqcv(W,s,c0,q,cv)
				rp(l)=getrp(getOmfromW(W,getcSig(c00,s0)),cv)
			enddo
		end subroutine

	end subroutine
	
	subroutine het
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		logical		:: iflag
		real		:: avrp
		real		:: rp(1)
		real		:: out(size(rp),nstates*ndraws)
		real		:: tab(3,size(rp),2)
		
		avc0=0.03
		do ilp=1,2
			lp=ilp-1
			print *,"ooooooooooooo het oooooooooooooooooooooooooooooooo"
			print *,"avc0=",avc0," lp=",lp
			out=-1; is=1
			do i=1,nstates
				do id=1,ndraws
					s=getstatelocs(i,lp,n)
					call getmrp
					
					out(:,is)=[rp]
					is=is+1
				enddo
			enddo
			call mdisp(out)
			do i=1,size(rp)
				tab(:,i,ilp)=quantile_v(out(i,:),[.05,.5,.95])
			enddo
			call mdisp(tab(:,:,ilp))
			call printtime
		enddo
		
	contains
		
		subroutine getmrp
			real, allocatable:: Om(:,:), W(:,:), Sig(:,:), Sig0(:,:)
			integer		:: i,j,l1,l2,ii
			real	:: maxs

			c0=getavc_c(avc0,0,s)			
			Sig=getcSig(c0,s)
			W=getWeigen(Sig,qmax+1)
			call setSCPCqcv(W,s,c0,q,cv)
			rp=-1
			do ii=1,2
				if(ii==1) then
					Sig0=getcSig(c0,s)
				else
					Sig0=eye(n)
				endif
				do i=1,2
					maxs=maxval(s(i,:),dim=1)
					do j=-1,1,2
						do l1=1,n
							do l2=1,n
								Sig(l1,l2)=Sig0(l1,l2)*exp(log(3.0)*j*(s(i,l1)+s(i,l2))/maxs)
							enddo
						enddo
						Om=matmul(transpose(W),matmul(Sig,W))
						rp=max(rp,getrp(Om,cv))
					enddo
				enddo
			enddo
		end subroutine

	end subroutine
	

	
	subroutine SunKim_distortion
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		logical		:: iflag
		real		:: avrp
		real		:: rp(1)
		real		:: out(1,nstates),gdens(n)
		real		:: tab(3,size(rp))
		out=-1; is=1
		do i=1,nstates
			do id=1,1
				s=getstatelocs(i,1.0,n,gdens)
				call getmrp
				out(:,is)=[rp]
				is=is+1
			enddo
		enddo
		print *,"oooooooooooooooooo SunKim_distortion ooooooooooooooooo"
		do i=1,size(rp)
			tab(:,i)=quantile_v(out(i,:),[.05,.5,.95])
		enddo
		
		call mdisp(out)
		
		call mdisp(tab)
		call printtime
		
	contains
		
		subroutine getmrp
			real, allocatable:: Om(:,:), W(:,:), Sig(:,:),Wg(:,:)
			integer, parameter	:: k1=1,k2=2,k=k1*k2+k1+k2
			integer			:: j

			W=getWSun(s,k1,k2)
			cv=getStudentcv(2*k,level)
			Wg=W
			do j=1,size(W,2)
				Wg(:,j)=Wg(:,j)*gdens
			enddo
			Om=matmul(transpose(Wg),W)
			rp=getrp(Om,cv)
		end subroutine

	end subroutine
	
	subroutine efficiency
		use plugin
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		real, allocatable	:: W(:,:)
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		logical		:: iflag
		real		:: avrp
		integer, parameter	:: nalts=3
		real		:: c1s(nalts),Sigs(n,n,nalts)
		real		:: rp(9,nalts)
		real		:: out(size(rp),nstates*ndraws)
		real		:: tab(3,size(rp),2)
		
		avc0=0.03
		do ilp=1,2
			lp=ilp-1
			print *,"ooooooooooooooooo efficiency oooooooooooooooooooooooooooo"
			print *,"avc0=",avc0," lp=",lp
			out=-1; is=1
			do i=1,nstates
				do id=1,ndraws
					s=getstatelocs(i,lp,n)
					call getmrp
					
					out(:,is)=[rp]
					is=is+1
				enddo
			enddo
			call mdisp(out)
			do i=1,size(rp)
				tab(:,i,ilp)=quantile_v(out(i,:),[.05,.5,.95])
			enddo
			call mdisp(tab(:,:,ilp))
			call savemat(trim(dir)//"efficiency"//convtos(nint(lp))//".txt",tab(:,:,ilp))
			call printtime
			
		enddo
		
	contains
	
		subroutine setlengths
			integer	:: ia
			do ia=1,nalts
				rp(ii,ia)=getlength(getOmfromW(W,Sigs(:,:,ia)),cv)
			enddo
			ii=ii+1
		end subroutine
		
		subroutine getmrp
			real, allocatable:: Sig(:,:), Walt(:,:), Wall(:,:),Om(:,:)
			real :: c1,lam0,cvalt
			integer	::qalt,i,ia
			
			c0=getavc_c(avc0,0,s)
			do ia=1,nalts
				if(ia==1) then
					Sigs(:,:,ia)=eye(n)
					c1s(ia)=-1
				else
					c1=merge(2,5,ia==3)*c0
					Sigs(:,:,ia)=getcSig(c1,s)
					c1s(ia)=c1
				endif
			enddo
			Sig=getcSig(c0,s)
			Wall=getWeigen(Sig,100)
			ii=1
			W=Wall(:,1:qmax+1)
			call setSCPCqcv(W,s,c0,q,cv)
			call setlengths
			
			W=Wall(:,1:qmax+1)
			call setSCPCqcv_alt(W,s,getOmfromW(W,Sigs(:,:,3)),c0,q,cv)
			call setlengths

			W=getWSun(s,1,2)
			cv=getcvfromOms(getOms(W,s,c0),size(W,2)-1,level)
			call setlengths
			
			W=getWKV(s,.3*getsmax(s))
			cv=getcvfromOms(getOms(W,s,c0),size(W,2)-1,level)
			call setlengths

			W=getclusterW(s,9)
			cv=getcvfromOms(getOms(W,s,c0),size(W,2)-1,level)
			call setlengths
			

			lam0=3*c0
			do ia=1,nalts
				rp(ii:ii+1,ia)=s2opt(s,c0,c1s(ia),lam0)
				Om=getOmfromW(Wall,Sigs(:,:,ia))
				call linds(Om(2:,2:),Om(2:,2:))
				rp(ii+2,ia)=sum(Om(2:,1)*matmul(Om(2:,2:),Om(2:,1)))/Om(1,1)
			enddo
			rp(1:ii,:)=rp(1:ii,:)/(2*gausscdfinv(1-level/2))
			rp(ii+3,:)=0.0/0.0
		end subroutine
	end subroutine

	subroutine tsrob
		real		::avclist(1)=[0.03]
		integer		:: nlist(5)=[50,100,200,500,1000]
		real, allocatable	:: s(:,:)
		integer		:: in,iavc,j,q
		real, allocatable:: Om(:,:), W(:,:), Sig(:,:)
		real	:: avc0,cv,c0
		
		print *,"oooooooo tsrob ooooooooooooooooooo"
		do in=1,size(nlist)
			n=nlist(in)
			if(allocated(s)) deallocate(s)
			allocate(s(1,n))
			do j=1,n
				s(:,j)=real(j)/n
			enddo
			do iavc=1,size(avclist)
				avc0=avclist(iavc)
				c0=getavc_c(avc0,0,s)			
				Sig=getcSig(c0,s)
				W=getWeigen(Sig,qmax+1)
				call setSCPCqcv(W,s,c0,q,cv)
				W(:,2:)=W(:,2:)*cv
				print *,n,iavc,q,getrobust_ts(W,c0,s)
			enddo
		enddo
	end subroutine
			
	subroutine MLErp
		use plugin
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		logical		:: iflag
		real		:: avrp
		real		:: rp(1)
		real		:: out(size(rp),nstates*ndraws)
		real		:: tab(3,size(rp),2)
		
		avc0=0.03
		do ilp=1,1
			lp=ilp-1
			print *,"ooooooooooooooooo MLErp oooooooooooooooooooooooooooo"
			print *,"avc0=",avc0," lp=",lp
			out=-1; is=1
			do i=1,nstates
				do id=1,ndraws
					s=getstatelocs(i,lp,n)
					call getmrp
					
					out(:,is)=[rp]
					is=is+1
				enddo
			enddo
			call mdisp(out)
			do i=1,size(rp)
				tab(:,i,ilp)=quantile_v(out(i,:),[.05,.5,.95])
			enddo
			call mdisp(tab(:,:,ilp))
			call printtime
		enddo
		
	contains
		
		subroutine getmrp
			real, allocatable:: Om(:,:), W(:,:), Sig(:,:)
			real	:: cv
			integer	:: i
			rp=getMLErp(s,avc0)
		end subroutine

	end subroutine

	
	subroutine newrob
		use optmod
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		real		:: rp(1)
		real		:: out(size(rp),nstates*ndraws)
		real		:: tab(3,size(rp),2)
		
		avc0=0.03
		do ilp=1,2
			lp=ilp-1
			print *,"ooooooooooooooooooooooooooooooooooooooooooooo"
			print *,"avc0=",avc0," lp=",lp
			out=-1; is=1
			do i=1,nstates
				do id=1,ndraws
					s=getstatelocs(i,lp,n)
					call getmrp
					
					out(:,is)=[rp]
					is=is+1
				enddo
			enddo
			call mdisp(out)
			do i=1,size(rp)
				tab(:,i,ilp)=quantile_v(out(i,:),[.05,.5,.95])
			enddo
			call mdisp(tab(:,:,ilp))
			call printtime
		enddo
		
	contains
		
		subroutine getmrp
			real, allocatable:: Om(:,:), W(:,:), Sig(:,:)
			real	:: cv
			c0=getavc_c(avc0,0,s)			
			Sig=getcSig(c0,s)
			W=getWeigen(Sig,qmax+1)
			call setSCPCqcv(W,s,c0,q,cv)		
			W(:,2:)=W(:,2:)*cv
			rp=getrobust_step(W,c0,s)
		end subroutine

	end subroutine
	
	subroutine SCPC_avc
		real		:: s(2,n),s0(2,n)
		real		:: c00,c0,cv
		integer		:: i,q,ilp,iavc,ii,is,p,l,id
		real		:: avc0,lp
		logical		:: iflag
		real		:: avrp
		real		:: rp(4)
		real		:: out(size(rp),nstates*ndraws)
		real		::avclist(4)=[0.003,0.01,0.03,0.10]
		real		:: tab(size(rp),3*size(avclist),2)
		tab=0
		do iavc=1,size(avclist)
			avc0=avclist(iavc)
			do ilp=1,2
				lp=ilp-1
				print *,"ooooooooooooooooooooooooooooooooooooooooooooo"
				print *,"avc0=",avc0," lp=",lp
				out=-1; is=1
				call rnset(123)
				do i=1,nstates
					do id=1,ndraws
						s=getstatelocs(i,lp,n)
						call getmrp
					
						out(:,is)=[rp]
						is=is+1
					enddo
				enddo
				do i=1,size(rp)
					tab(i,3*(iavc-1)+1:3*iavc,ilp)=quantile_v(out(i,:),[.05,.5,.95])
				enddo
				call mdisp(tab(:,:,ilp))
				call savemat(trim(dir)//"SCPCtab"//convtos(nint(lp))//".txt",tab(:,:,ilp))
				
				call printtime
			enddo
		enddo		
	contains
		
		subroutine getmrp
			real, allocatable:: Om(:,:), W(:,:), Sig(:,:)
			real	:: cv
			c0=getavc_c(avc0,0,s)			
			Sig=getcSig(c0,s)
			W=getWeigen(Sig,merge(55,qmax,avc0<0.01))
			call setSCPCqcv(W,s,c0,q,cv)
			rp(1)=q
			rp(2)=100*(log(2.0)/c0)/getsmax(s)
			rp(3)=getIlength(q,cv)/(2*gausscdfinv(1-level/2))
			rp(4)=getrp(getOmfromW(W,Sig),cv)
		end subroutine

	end subroutine


end module
	
program mfort
	use runmod
	implicit none
	integer	:: i
	
	call rnset(15)
!$omp parallel do
	do i=1,100
		call rnset(3491084*i)
	enddo
	call inittime
	
	call mkGQxw(GQxw)
	
	call tsrob
	n=1000
	call SunKim_distortion
	n=500
	call efficiency
	call newrob
	call SCPC_avc
	call het
	call locerror
	call SunKim_distortion		
	call MLErp
end program