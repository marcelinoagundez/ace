
c      user subroutine for mnewt
c      it evaluates the functions F_i(x_j)=0 and the jacobian

       subroutine usrfun(x,n,NP,fvec,fjac)
       
       include 'ace.common'       
       integer i,j,k,l,NP,n
       real*8 x(NP),fvec(NP),fjac(NP,NP),s1,s2,s3,s4,s5

       if(n.ne.neqt)stop ' E- Error in mnewt usrfun'

c      ...get current values of pi lagrange multipliers
       do i=1,npilag
          pilag(i)=x(i)
       enddo
       do l=1,ncondensed
          dnc(l)=x(npilag+l)
       enddo
       dlnn=x(neqt)

c      ...evaluate functions f_i(x_j)
       do k=1,npilag                             ! element+charge conservation equations
          s1=0.0d0
          s2=0.0d0
          s3=0.0d0
          s4=0.0d0
          s5=0.0d0
          do j=1,nspec
             s4=s4+nat(j,k)*abun(j)
             if(.not.condensed(j))then
                do i=1,npilag
                   s1=s1+nat(j,k)*nat(j,i)*abun(j)*pilag(i)
                enddo
                s3=s3+nat(j,k)*abun(j)
                s5=s5+nat(j,k)*abun(j)*mu(j)
             else
                s2=s2+nat(j,k)*dnc(j-ngas)
             endif
          enddo
          fvec(k)=s1+s2+s3*dlnn-b0(k)+s4-s5
       enddo
       do l=1,ncondensed                         ! condensed species equations
          j=ngas+l
          s1=0.0d0
          do i=1,npilag
             s1=s1+nat(j,i)*pilag(i)
          enddo
          fvec(npilag+l)=s1-mu(j)
       enddo
       s1=0.0d0                                  ! total mole conservation equation
       s2=0.0d0
       s3=0.0d0
       do j=1,ngas
          do i=1,npilag
             s1=s1+nat(j,i)*abun(j)*pilag(i)
          enddo
          s2=s2+abun(j)
          s3=s3+abun(j)*mu(j)
       enddo
       fvec(neqt)=s1+(s2-abuntot)*dlnn-abuntot+s2-s3

c      ...evaluate derivatives of functions d[f_i(x_j)]/dx_j
       do k=1,npilag                             ! f_i(x_j) = element+charge
          do i=1,npilag
             fjac(k,i)=0.0d0
             do j=1,nspec
                fjac(k,i)=fjac(k,i)+nat(j,k)*nat(j,i)*abun(j)
             enddo
          enddo
          do l=1,ncondensed
             j=ngas+l
             fjac(k,npilag+l)=nat(j,k)
          enddo
          fjac(k,neqt)=0.0d0
          do j=1,ngas
             fjac(k,neqt)=fjac(k,neqt)+nat(j,k)*abun(j)
          enddo
       enddo
       do l=1,ncondensed                         ! f_i(x_j) = condensed
          j=ngas+l
          do i=1,npilag
             fjac(npilag+l,i)=nat(j,i)
          enddo
          do i=npilag+1,neqt
             fjac(npilag+l,i)=0.0d0
          enddo
       enddo
          do i=1,npilag                          ! f_i(x_j) = total mole
             fjac(neqt,i)=0.0d0
             do j=1,ngas
                fjac(neqt,i)=fjac(neqt,i)+nat(j,i)*abun(j)
             enddo
          enddo
          do l=1,ncondensed
             fjac(neqt,npilag+l)=0.0d0
          enddo
          fjac(neqt,neqt)=0.0d0
          do j=1,ngas
             fjac(neqt,neqt)=fjac(neqt,neqt)+abun(j)
          enddo
          fjac(neqt,neqt)=fjac(neqt,neqt)-abuntot

       return
       end
