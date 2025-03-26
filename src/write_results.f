
       subroutine write_results
       
       include 'ace.common'       
       integer i,j,out,ncm,i0,iz,ne,ilen
       real*8 wabun(nmaxspec)
       character*(nmaxcharfile) outfile,zabfile
       character*40 txtatom

c     ...if single point calculation, write .spec input file for 0d
       if(napt.eq.1)then
          if(init)return
          outfile(1:nmaxcharfile)=''
          outfile(1:ncharmodel)=namemodel(1:ncharmodel)
          outfile(ncharmodel+1:ncharmodel+9)='_out.spec'
          open(unit=9,err=110,file=outfile,status='unknown')
          write(9,1000)specfile,thermfile,aptfile
          do i=1,nspec
             wabun(i)=max(abun(i)/abuntot,abunmin)    ! evaluate molar fraction
             ne=0
             txtatom(1:40)=''
             ilen=0
             do j=1,nelem                             ! build up atoms text
                if(nat(i,j).eq.0)cycle
                ne=ne+1
                txtatom(ilen+1:ilen+2)=elem(j)(1:2)
                ilen=ilen+2
                if(elem(j)(2:2).eq.' ')ilen=ilen-1
                txtatom(ilen+1:ilen+2)='  '
                ilen=ilen+2
                if(nat(i,j).lt.10)then
                   write(txtatom(ilen+1:ilen+1),'(i1)')nat(i,j)
                   ilen=ilen+1
                else
                   write(txtatom(ilen+1:ilen+2),'(i2)')nat(i,j)
                   ilen=ilen+2
                endif
                txtatom(ilen+1:ilen+3)='   '
                ilen=ilen+3
             enddo
             ilen=ilen-3
             write(9,1010)spec(i),wabun(i),-99999.0e0,-999.0e0,0.0e0,
     #       0.0e0,charge(i),nattot(i),ne,txtatom
          enddo
          close(9)
          return
       endif

c     ...initialize output files
       if(init)then
          noutfile=nspec/nspecxfile
          if(noutfile*nspecxfile.lt.nspec)noutfile=noutfile+1
          if(noutfile.gt.12)goto 100
          out=9
          ncm=ncharmodel
          do i=1,noutfile
             outfile(1:nmaxcharfile)=''
             outfile(1:ncm)=namemodel(1:ncm)
             outfile(ncm+1:ncm+11)='_out_00.dat'
             if(i.lt.10)write(outfile(ncm+7:ncm+7),'(i1)')i
             if(i.ge.10)write(outfile(ncm+6:ncm+7),'(i2)')i
             open(unit=out,err=110,file=outfile,status='unknown')
             write(out,1005)specfile,thermfile,aptfile
             i0=nspecxfile*(i-1)+1
             iz=min0(nspecxfile*i,nspec)
             write(out,1020)(j,j=1,iz-i0+1+3)
             write(out,1030)(spec(j),j=i0,iz)
             out=out+1
          enddo
          zabfile(1:nmaxcharfile)=''
          zabfile(1:ncm)=namemodel(1:ncm)
          zabfile(ncm+1:ncm+4)='.zab'
          open(unit=8,err=120,file=zabfile,status='unknown')
          write(8,1025)napt,nspec
          write(8,1035)(spec(j),j=1,nspec)

c     ...write output abundances
       else
          do i=1,nspec
             wabun(i)=dmax1(abun(i)/abuntot,abunmin)  ! evaluate molar fraction
          enddo
          out=9
          do i=1,noutfile
             i0=nspecxfile*(i-1)+1
             iz=min0(nspecxfile*i,nspec)
             write(out,1040)altitude,pressure,tk,
     #       (dlog10(wabun(j)),j=i0,iz)
             out=out+1
          enddo
          write(8,1040)altitude,pressure,tk,
     #    (dlog10(wabun(j)),j=1,nspec)
       endif

c     ...close output files
       if(iapt.eq.napt)then
          out=9
          do i=1,noutfile
             close(out)
             out=out+1
          enddo
          close(8)
       endif

       return

100    write(*,*)' E- Number of output files > 12'
       stop
110    write(*,*)' E- Error opening outfile ',outfile
       stop
120    write(*,*)' E- Error opening .zab outfile ',zabfile
       stop

1000   format(
     # '!',/,
     # '! thermochemical equilibrium model:',/,
     # '! species file : ',a,/,
     # '! therm file :   ',a,/,
     # '! (p,T) file :   ',a,/,
     # '! ..............................................................
     #...................................',/,
     # '! Abundances expressed as decimal logarithm of molar fraction:',
     # /,'!',/,
     # '!species      molar fr    DHf:kJ/mol   mu:Debye   polar:A^3',
     # '    IP:eV  +/- Natom Nelem    atoms',/,
     # '!------------ ---------------------------------------------',
     # '--------------------------    ','--------------------------',
     # '---------------')
1005   format(
     # '!',/,
     # '! thermochemical equilibrium model:',/,
     # '! species file : ',a,/,
     # '! therm file :   ',a,/,
     # '! (p,T) file :   ',a,/,
     # '! ..............................................................
     #...................................',/,
     # '! Abundances expressed as decimal logarithm of molar fraction:'
     # )
1010   format(1x,a15,1x,1pg9.3,2x,0pf10.3,2x,0pf10.3,1x,0pf10.6,
     # 0pf10.6,3x,i2,4x,i2,4x,i2,4x,a)
1020   format('!',5x,i2,12x,i2,12x,i2,5x,9999(4x,i4,8x))
1025   format(1x,i4,36x,'! Number of heights',/,
     #        1x,i4,36x,'! Number of species')
1030   format('    location    pressure[bar]    Tk[K]    ',
     # 2x,9999(1x,a15))
1035   format('    location    pressure[bar]    Tk[K]    ',
     # 2x,9999(a15,1x))
1040   format(1x,1pg12.6,2x,1pe12.4,2x,0pf10.2,9999(2x,0pf12.6,2x))

       end
