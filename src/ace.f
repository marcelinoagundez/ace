
       include 'ace.common'
       integer i

       write(*,1000)
       open(unit=1,err=100,file='ace.inp',status='old')

       i=0
       do
          call reset_variables
          read(1,*,err=110,end=10)inpfile
          if(inpfile(1:3).eq.'END'.or.inpfile(1:3).eq.'end')exit
          i=i+1
          if(i.gt.nmaxmodel)stop' E- Too many models'
          call read_input_files
          call read_spec_therm
          call compute_chemical_equilibrium
       enddo
10     close(1)

       write(*,1010)
       stop

100    write(*,*)' E- Error opening ace.inp'
       stop
110    write(*,*)' E- Error reading ace.inp'
       stop

1000   format(/,
     # ' **********************************************************',/,
     # ' *                           ACE                          *',/,
     # ' *        Atmospheric chemical Equilibrium program        *',/,
     # ' *                                                        *',/,
     # ' *                    version May 2023                    *',/,
     # ' **********************************************************',/)
1010   format(/,
     # ' ==========================================================',/,
     # ' =                    Calculations done                   =',/,
     # ' =                       Exiting ACE                      =',/,
     # ' ==========================================================',/)

       end




c_______________________________________________________________________

       subroutine reset_variables
       
       include 'ace.common'       

c      Input/output files
       inpfile(1:nmaxcharfile)=''
       namemodel(1:nmaxcharfile)=''
       ncharmodel=0
       specfile(1:nmaxcharfile)=''
       thermfile(1:nmaxcharfile)=''
       aptfile(1:nmaxcharfile)=''
       noutfile=0

c      Physical parameters
       napt=0
       a_apt(:)=0.0d0
       p_apt(:)=0.0d0
       t_apt(:)=0.0d0
       altitude=0.0d0
       pressure=0.0d0
       tk=0.0d0

c      Elements data
       nelem=0
       elem(:)(1:2)=''
       ilenelem(:)=0
       elfab(:)=0.0d0
       mat(:)=0.0d0
       zat(:)=0

c      Species data
       nspec=0
       spec(:)(1:nmaxcharspec)=''
       ilenspec(:)=0
       nattot(:)=0
       nat(:,:)=0
       idzat(:,:)=0
       charge(:)=0
       condensed(:)=.false.
       abun(:)=0.0d0
       abunmax(:)=0.0d0
       mu(:)=0.0d0
       abuntot=0.0d0
       id_electr=0
       ngas=0
       ncondensed=0

c      Therm data
       type_therm(:)=0
       natherm(:)=0
       ntemp_therm(:)=0
       temp_therm(:,:)=0.0d0
       atherm(:,:,:)=0.0d0
       spec_tkout(:)=.false.
       spec_warning(:)=.false.

c      Numerical and model parameters
       neqt=0
       npilag=0
       b0(:)=0.0d0
       b0max=0.0d0
       b0min=0.0d0
       pilag(:)=0.0d0
       dlnn=0.0d0
       dnc(:)=0.0d0
       ion=.false.
       converge=.false.

       return
       end
