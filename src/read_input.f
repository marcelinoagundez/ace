
       subroutine read_input_files
       
       include 'ace.common'       
       integer i
       character*100 line

       i=nmaxcharfile
       do while (inpfile(i:i).ne.'.')
          i=i-1
       enddo
       ncharmodel=i-1
       namemodel(1:ncharmodel)=inpfile(1:ncharmodel)

       write(*,1000)inpfile
       open(unit=2,err=100,file=inpfile,status='old')
       read(2,*)aptfile
       open(unit=3,err=110,file=aptfile,status='old')
       napt=0
       do i=1,nmaxapt
          line(1:100)=''
          read(3,'(a100)',end=10)line
          if(line(1:1).eq.'!')cycle
          napt=napt+1
          read(line(1:100),*,err=120)
     #    a_apt(napt),p_apt(napt),t_apt(napt)
       enddo
10     close(3)
       read(2,*)specfile
       read(2,*)thermfile
       close(2)

       return

100    write(*,*)' E- Error opening model input file: ',inpfile
       stop
110    write(*,*)' E- Error opening (p,T) file: ',aptfile
       stop
120    write(*,*)' E- Error reading, in (p,T) file, line:'
       write(*,*)line
       stop

1000   format(2x,'I- Running model ',a)

       end




c_______________________________________________________________________

       subroutine read_spec_therm
       
       use ace_data
       include 'ace.common'       
       integer i,j,k,l,ntherm,telfab,
c      telfab:   type of input elemental abundances
c                =1: relative abundance
c                =2: logarithmic abundance relative to H=12
     # read_nattot,read_nelem,read_charge,read_type,
     # read_nat(nmaxelemxs),read_ntemp,read_na
       real*8 read_a(9,nmaxtemp_therm),read_temp(nmaxtemp_therm+1)
       character*256 line
       character*50 read_spec
       character*2 read_txt(2*nmaxelemxs),read_elem(nmaxelemxs)
       logical fnd_spec(nmaxspec),fnd_read_elem,idem

c     ...open species file
       write(*,1000)specfile
       open(unit=3,err=100,file=specfile,status='old')

c     ...read number of elements
       do 
          line(1:256)=''
          read(3,'(a256)',end=10)line
          if(line(1:1).eq.'!')cycle
          read(line,*,err=110)nelem,telfab
          exit
       enddo
       if(nelem.gt.nmaxelem)then
          write(*,*)' E- Number of elements too high'
          stop
       endif

c     ...read elemental abundances
       i=0
       do while(i.lt.nelem)
          line(1:256)=''
          read(3,'(a256)',end=10)line
          if(line(1:1).eq.'!')cycle
          i=i+1
          read(line,*,err=110)elem(i),elfab(i)
          ilenelem(i)=2
          if(elem(i)(2:2).eq.' ')ilenelem(i)=1
          do j=1,nmaxelem
             call same_string(2,elem(i),element(j),idem)
             if(idem)then
                mat(i)=atomic_weight(j)
                zat(i)=j
             endif
          enddo
          if(mat(i).eq.0.0d0)then
             write(*,*)' E- Element ',elem(i)(1:2),' not recognized'
             stop
          endif
       enddo
       if(telfab.eq.2)then                 ! convert logarithmic to relative abundance
          do i=1,nelem
             elfab(i)=10.0d0**(elfab(i)-12.0d0)
          enddo
       endif

c     ...read species
       nspec=0
       ncondensed=0
       do while(nspec.le.nmaxspec+1)
          line(1:256)=''
          read(3,'(a256)',end=10)line
          if(line(1:1).eq.'!')cycle
          nspec=nspec+1          
          if(nspec.gt.nmaxspec) stop ' E- Too many species'
          read_spec(1:50)=''
          read(line,*,err=110,end=10)read_spec
          l=50
          do while (read_spec(l:l).eq.' ')
             l=l-1
          enddo
          spec(nspec)(1:l)=read_spec(1:l)
          ilenspec(nspec)=l
          if(ilenspec(nspec).gt.nmaxcharspec)then 
             write(*,*)' E- Many characters for: ',spec(nspec)
             stop
          endif
          call same_string(3,spec(nspec)(l-2:l),'(s)',idem) ! solid species
          if(idem)condensed(nspec)=.true.
          call same_string(3,spec(nspec)(l-2:l),'(l)',idem) ! liquid species
          if(idem)condensed(nspec)=.true.
          call same_string(3,spec(nspec)(l-2:l),'(c)',idem) ! generic condensed species
          if(idem)condensed(nspec)=.true.
          if(condensed(nspec))ncondensed=ncondensed+1
          if(ncondensed.gt.0.and..not.condensed(nspec))
     #    stop ' E- Condensed species must be at end of file .spec'
       enddo
10     close(3)
       ngas=nspec-ncondensed
       if(nspec+ncondensed.gt.nmaxspec)stop' E- Need larger nmaxspec'

c     ...read therm data
       write(*,1010)thermfile
       open(unit=3,err=120,file=thermfile,status='old')
       ntherm=0
       do i=1,nspec
          fnd_spec(i)=.false.
       enddo
       do while(ntherm.le.nmaxspec)

          line(1:256)=''                   ! read parameters from therm file
          read(3,'(a256)',end=20)line
          if(line(1:1).eq.'!')cycle
          read(line,*,err=130,end=20)read_spec,read_nattot,read_nelem,
     #    read_charge,read_type
          read(3,*,err=135)(read_txt(j),j=1,2*read_nelem)
          do i=1,2*read_nelem
             if(((-1)**i).lt.0)read_elem((i+1)/2)(1:2)=read_txt(i)(1:2)
             if(((-1)**i).gt.0)read(read_txt(i)(1:2),*)read_nat(i/2)
          enddo
          read(3,*,err=140)read_ntemp,(read_temp(j),j=1,read_ntemp+1)
          if(read_type.eq.1)read_na=7
          if(read_type.eq.7)read_na=7
          if(read_type.eq.9)read_na=9
          if(read_type.ne.1.and.read_type.ne.7.and.read_type.ne.9)then
             write(*,*)' E- unknown therm data type for: ',read_spec
             stop
          endif
          do i=1,read_ntemp
             read(3,*,err=140)(read_a(j,i),j=1,read_na)          
          enddo

          i=50
          do while (read_spec(i:i).eq.' ')
             i=i-1
          enddo
          if(i.gt.nmaxcharspec)then 
             write(*,*)' E- Many characters for (therm): ',read_spec
             stop
          endif

          do i=1,nspec                      ! verify whether species is included or not
             if(fnd_spec(i))cycle           ! ... and if so, save parameters:
             call same_string(nmaxcharspec, ! - number of atoms: nattot(i)
     #       read_spec,spec(i),idem)        ! - number of atoms of each element: nat(i,j)
             if(idem)then                   ! - charge: charge(i)
                ntherm=ntherm+1             ! - thermo data:    - type_therm(i)
                nattot(i)=read_nattot                         ! - ntemp_therm(i)
                if(read_nelem.gt.nmaxelemxs-1)then            ! - temp_therm(i,k)
                   write(*,*)' E- Many elements in: ',spec(i) ! - atherm(i,j,k)
                   stop
                endif
                do j=1,read_nelem
                   fnd_read_elem=.false.
                   do k=1,nelem
                      call same_string(2,read_elem(j),elem(k),idem)
                      if(idem)then
                         nat(i,k)=read_nat(j)
                         idzat(i,k)=zat(k)
                         fnd_read_elem=.true.
                         exit
                      endif
                   enddo                   
                   if(.not.fnd_read_elem)then
                      write(*,*)' E- Error on elements for: ',spec(i)
                      stop
                   endif
                enddo
                charge(i)=read_charge
                type_therm(i)=read_type
                natherm(i)=read_na
                ntemp_therm(i)=read_ntemp
                do k=1,ntemp_therm(i)+1
                   temp_therm(i,k)=read_temp(k)
                   if(k.eq.ntemp_therm(i)+1)exit
                   do j=1,natherm(i)
                      atherm(i,j,k)=read_a(j,k)
                   enddo
                enddo
                fnd_spec(i)=.true.
                exit
             endif
          enddo

          read_spec(1:50)=''                ! reset 'read_...' variables
          read_nattot=0
          read_nelem=0
          read_charge=0
          read_type=0
          read_na=0
          read_txt(:)(1:2)=''
          read_elem(:)(1:2)=''
          read_nat(:)=0
          read_ntemp=0
          read_temp(:)=0.0d0
          read_a(:,:)=0.0d0

       enddo
20     close(3)

c     ...verify that all species have thermo data
       if(ntherm.ne.nspec)then
          write(*,*)' E- Error searching for thermo data'
          stop
       endif

c     ...verify that the total number of atoms of each species is correct
       do i=1,nspec
          k=0
          do j=1,nelem
             k=k+nat(i,j)
          enddo
          if(k.ne.nattot(i))then
             write(*,*)' E- Wrong number of atoms for: ',spec(i)
             stop
          endif
       enddo

c     ...verify that a species does not appear more than once
       do i=1,nspec-1
          l=ilenspec(i)
          do j=i+1,nspec
             k=ilenspec(j)
             if(l.ne.k)cycle
             if(spec(i)(1:l).eq.spec(j)(1:l))then
                write(*,*)' E- Twice species ',spec(i)(1:l)
                stop
             endif
          enddo
       enddo

c     ...verify species charge
       j=0
       k=0
       ion=.false.
       do i=1,nspec
          if(charge(i).ne.0)ion=.true.
          if(charge(i).gt.0)j=j+1
          if(charge(i).lt.0)k=k+1
          if(charge(i).ne.0.and.condensed(i))then
             write(*,*)' E- Charge in condensed species ',spec(i)
             stop
          endif
       enddo
       if(ion.and.(j.eq.0.or.k.eq.0))then
          write(*,*)' E- Missing (+/-) charged species'
          stop
       endif

c     ...identify electron species if present
       if(ion)then
          id_electr=0
          do i=1,nspec
             if(nattot(i).eq.0.and.charge(i).eq.-1)id_electr=i
          enddo
       endif

c     ...verify that number of atoms is non zero except for electron
       do i=1,nspec
          if(ion.and.i.eq.id_electr)cycle
          if(nattot(i).eq.0)then
             write(*,*)' E- Number of atoms zero for: ',spec(i)
             stop
          endif
       enddo

       write(*,1020)nspec
       write(*,1030)ncondensed

       return

100    write(*,*)' E- Error opening species file'
       stop
110    write(*,*)' E- Error reading species file at line:'
       write(*,*)'  ',line
       stop
120    write(*,*)' E- Error opening therm file'
       stop
130    write(*,*)' E- Error reading therm file at line:'
       write(*,*)'  ',line
       stop
135    write(*,*)' E- Error reading therm file elements for: ',
     # read_spec
       stop
140    write(*,*)' E- Error reading thermo data for: ',read_spec
       stop

1000   format(2x,'I- Reading species from ',a45)
1010   format(2x,'I- Reading therm data from ',a45)
1020   format(2x,'I- Data read for ',i4,' species')
1030   format(2x,'I- ',i4,' condensed species')

       end




c_______________________________________________________________________

       subroutine same_string(nmaxchar,string1,string2,idem)
c      ...evaluate if two character strings are similar (case-sensitive) 
c            idem=.false. --> different strings
c            idem=.true.  --> similar strings

       implicit none

       integer nmaxchar,l1,l2,k,l
       character*(nmaxchar) string1,string2,cwk1,cwk2
       logical idem

       idem=.false.

       l1=nmaxchar
       do while(string1(l1:l1).eq.' ')
          l1=l1-1
       enddo
       l2=nmaxchar
       do while(string2(l2:l2).eq.' ')
          l2=l2-1
       enddo
       if(l1.eq.l2)then
          cwk1(1:nmaxchar)=''
          cwk2(1:nmaxchar)=''
          do k=1,l1
             l=ichar(string1(k:k))
             cwk1(k:k)=char(l)
             l=ichar(string2(k:k))
             cwk2(k:k)=char(l)
          enddo
          if(cwk1(1:l1).eq.cwk2(1:l1))idem=.true.
       endif

       return
       end
