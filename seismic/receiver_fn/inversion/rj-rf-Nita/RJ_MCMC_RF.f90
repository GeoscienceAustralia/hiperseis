!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This program does a Transdimensional inversion of Receiver functions
! with the reversible jump algorithm

! Thomas Bodin, ANU, December 2011
!*********************************************************************
!/////////////////////////////////////////////////////////////////////
!********************************************************************

program RJ_MCMC_RF
use mpi
implicit none

!             ----------------------------------
!             BEGINING OF THE USER-DEFINED PARAMETERS
!             ----------------------------------


!-----------------------------------------
! Parameters of the Markov chain
!-----------------------------------------

integer, parameter :: burn_in = 650000  !Burn-in period
integer, parameter :: nsample = 950000  !Post burn-in
integer, parameter :: thin = 100        !Thinning of the chain

! Each chain is run for 'burn_in + nsample' steps in total. The first burn-in samples are discarded as burn-in steps,
! only after which the sampling algorithm is assumed to have converged. To eliminate dependent samples in the ensemble
! solution, every thinn model visited in the second part is selected for the ensemble.
! 
! The convergence of the algorithm is monitored with a number of indicators such as acceptance rates, and sampling
! efficiency is optimized by adjusting the variance of the Gaussian proposal functions 


!------------------------------------------------
! Prior distribution (Bounds odf the model space)
!------------------------------------------------

integer,parameter  :: npt_min = 2 !minimun number of layers is 2. 1 does not work.
integer, parameter :: npt_max = 20

!depth
real, parameter :: d_min = 0
real, parameter :: d_max = 60

real, parameter :: beta_min = 3.5 ! Mean of the Uniform prior on Vs
real,parameter  :: beta_max = 3.5 ! Mean of the Uniform prior on Vs
real, parameter :: width = 2.0 ! Lower and upper bound of the prior are [mean-theta , mean+theta] Note that beta_min needs to be equal to beta_max.

double precision, parameter ::    Ar_max=0.1000  !Upper bound for noise parameter
double precision, parameter ::    Ar_min=0.00100 !Lower bound for noise parameter 


!-----------------------------------------
! Sdt for Proposal distributions
!-----------------------------------------

!The convergence of the algorithm is monitored with a number of indicators such as
! acceptance rates, and sampling efficiency is optimized by adjusting the variance of
! the Gaussian proposal functions 

! These values have to be "tuned" so the acceptance rates written in OUT/mpi.out
! are as close as possible to 44%. This determined the efficiency of posterior sampling.
!  If AR larger than 44%, increase the Sdt for less Acceptance.
!  If AR_* smaller than 44%, decrease the Sdt for more

real, parameter :: pd1 = 0.25     !proposal on change in position  
real, parameter :: pv1 = 0.05     !proposal on velocity
real, parameter :: pd2 = 0.75     !proposal on change in position 
real, parameter :: pv2 = 0.15     !proposal on velocity
real, parameter :: sigmav = 0.06  !proposal on velocity when Birth move  
real, parameter :: pAr = 0.002     !proposal for change in noise parameter 


!--------------------------------------------
! Parameters for Displaying results 
!-------------------------------------------- 

integer, parameter :: display = 20000 ! display results in OUT/mpi.out
!every display samples
double precision, parameter ::    sig=2.5 !fs/a

 !discretezation for the posterior distribution.
 !specifies the number of velocity pixels in (Vs, depth)
 !THis is becuas ethe solution is visualized as an histogram, 
 !and hence we need to define the number od bins

integer, parameter :: disd = 100 !depth
integer, parameter :: disv = 75  !velocity
integer, parameter :: disA = 75 !for noise parameter 

!depth of model for display
real, parameter :: prof = 60 


!--------------------------------------------
! Parameters for Receiver Function
!-------------------------------------------- 

!Those parameters are associated with the data vector which is a waveform

real, parameter :: fs = 6.25        !sampling frequency
real, parameter :: gauss_a = 2.5    !number 'a' defining the width of 
!the gaussian filter in the deconvolution  
real, parameter :: water_c = 0.0001 !water level for deconvolution 
real, parameter :: angle = 31.85       ! angle of the incoming wave 
real, parameter :: time_shift = 5   !time shift before the first p pulse
integer, parameter :: ndatar = 157  !Number of data points

!Since we need to compute numerically the inverse of the data covariance matrix numerically,
!please try to keep ndatar as small as possible (i.e. fs as small as possible), without losing
!information on the waveform of course.

integer, parameter :: maxdata = 250 !needs to be > than ndatar

real, parameter ::    pi = 3.14159265      !don't touch
real, parameter ::    rad = 0.017453292    !don't touch
real, parameter ::    v60 = 8.043      ! Vs at 60km depth 
!v60 is needed to compute the ray parameter from angle. 
! The equations used later in the code are :
!ppara=sin(angle*rad)/v60
!din=asin(ppara*beta(npt)*vpvs(npt))/rad


!--------------------------------------------
! Parameter for truncation of eigenvalues of Cd
!--------------------------------------------

! This will stabiolie the inversion of Cd
! Choose this number by looking at picard plot from Picard.m


double precision, parameter :: wmin=0.000001

!             ----------------------------------
!             END OF THE USER -DEFINED PARAMETERS
!             ----------------------------------



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!***************************************************************

! DECLARATIONS OF VARIABLES

!****************************************************************

real , EXTERNAL    ::    gasdev,ran3
real log, sqrt

integer i,npt,sample,ind,accept,l,th,ount,npt_prop,nb
integer birth, death,move,value,noisr,noisd,ind2,v,j,chwidd,chwidr
integer histo(npt_max),histos(npt_max),histoch(disd),histochs(disd)
real beta(npt_max),h(npt_max),vpvs(npt_max),qa(npt_max)
real qb(npt_max),av(disd),post(disd,disv),PrP(2),PrV(2),AcP(2),AcV(2),mo(disd)
real PrB,AcB,PrD,AcD,Prnr,Acnr,avs(disd),posts(disd,disv),out,LSr,LSr_prop,LSr_min
real ppara, din,ht,d,logprob,convAr(nsample+burn_in),convArs(nsample+burn_in)
real wdata(maxdata),best_datar(maxdata),d_obsr(maxdata),voro(npt_max,3)
real voro_prop(npt_max,3),conv(nsample+burn_in),ncell(nsample+burn_in)
real convs(nsample+burn_in),ncells(nsample+burn_in),t1,t2
real like,like_prop,liker,liker_prop,u,pv_min,pv_max,like_min,xi
real ML_Ar(disA),ML_Ars(disA)
double precision br(ndatar),Ar,Ar_prop,logrsig
double precision AIr(ndatar,ndatar),AIr_prop(ndatar,ndatar),Cd(ndatar,ndatar),UU(ndatar,ndatar),W(ndatar),x(ndatar)
double precision vv(ndatar,ndatar),bes(ndatar),b(ndatar)
double precision r

!for MPI
integer ra,ran,rank, nbproc, ierror

!Add by Sheng RSES ANU Aug-2018
integer clock, seed

! Variables for command line interface
integer :: num_cli_args
integer :: cli_status
integer :: dummy_cli_len
character(len=256) :: input_file
character(len=256) :: output_folder

!***********************************************************************

!**************************************************************

!                  CHECK AND READ ARGUMENTS

!**************************************************************
Acnr=0.
liker_prop=0
lsr_prop=0
Prnr=0

num_cli_args = command_argument_count()
if (num_cli_args /= 2) then
  print *, 'Usage: run input_filename output_folder'
  call exit()
end if

call get_command_argument(1, input_file, dummy_cli_len, cli_status)
if (cli_status < 0) then
  print *, 'Input file path too long'
  call exit()
else if (cli_status > 0) then
  print *, 'Unknown error retrieving input filename'
  call exit()
else
  print *, 'Input file: ', TRIM(input_file)
end if

call get_command_argument(2, output_folder, dummy_cli_len, cli_status)
if (cli_status < 0) then
  print *, 'Output path too long'
  call exit()
else if (cli_status > 0) then
  print *, 'Unknown error retrieving output path'
  call exit()
else
  print *, 'Output path: ', TRIM(output_folder)
end if


!**************************************************************

!                  INITIALIZATION

!**************************************************************

CALL cpu_time(t1)  !Tic. start counting time 
 

 !Start Parallelization of the code. From now on, the code is run on each
 !processor independently, with ra = the number of the proc.
 ! Note that running in parallel does not reduce total CPU time, instead
 ! it just samples more random states and averages more samples at the end.
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  write (*,*) nbproc, rank 

  nb=nbproc
  ra=rank
  ran=rank
 conv=0
 ncell=0

! Add by Sheng RSES ANU Aug-2018
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + rank
ra=int(seed)
write(*,*)'Rank',rank,'using seed',ra


!**************************************************************

!                  READ INPUT FILE

!**************************************************************


open(55, file=TRIM(input_file), status='old')
do i=1,ndatar
read(55,*)u,d_obsr(i)
end do
close(55)! close the file


!**************************************************************

!    Draw the first model randomly from the prior distribution

!**************************************************************


! Initilalise sigma (noise parameter)
Ar =  Ar_min+ran3(ra)*(Ar_max-Ar_min)!

! Initial number of cells
npt = npt_min+int(ran3(ra)*(8-npt_min))

! We place them randomly
do i=1,npt
	voro(i,1)=d_min+ran3(ra)*(d_max-d_min)
	call priorvalue(voro(i,1),d_min,d_max,beta_min,beta_max,width,pv_min,pv_max)
	voro(i,2)=pv_min+(pv_max-pv_min)*ran3(ra)
	voro(i,3)=1.73!vpvs_min+ran3(bb)*(vpvs_max-vpvs_min)
enddo

do i=npt+1,npt_max
voro(i,1)=0
voro(i,2)=0
voro(i,3)=1.73
enddo

! Add by Sheng RSES ANU Aug-2018
! add -DDEBUG in compilatoin to turn on following check
! #ifdef DEBUG
!     write(*, *) "Initial Model: num of layers: ", npt
! 	do i=1,npt
!         write(*, *) "voro(", i, ")=", voro(i,1), voro(i,2), voro(i,3)
! 	enddo
! #endif

call voro2qmodel(voro,npt,npt_max,d_min,d_max,beta,h,vpvs,qa,qb)

! Add by Sheng RSES ANU Aug-2018
! add -DDEBUG in compilatoin to turn on following check
! #ifdef DEBUG
!     write(*, *) "Voro Model: num of layers: ", npt
! 	do i=1,npt
!         write(*, *) i, "beta=", beta(i), "h=", h(i), "vpvs=", vpvs(i), "qa=", qa(i), "qb=", qb(i)
! 	enddo
!     write(*, *) "Start Iteration"
! #endif

! Compute Rf for initial model
ppara=sin(angle*rad)/v60
din=asin(ppara*beta(npt)*vpvs(npt))/rad
call theo(&
	npt, beta, h, vpvs, qa, qb, fs, din,&
	gauss_a, water_c, time_shift, ndatar, wdata )


!***************************************************

!       Get the data covariance Matrix Cd for RF

!***************************************************

r=exp(-1/(2*sig**2))
do i=1,ndatar
	do j=1,ndatar
		Cd(i,j)=r**abs(i-j)
	enddo		
enddo

!***************************************************

!       SVD Decomposition of the Matrix Cd

!***************************************************


UU=Cd
 call dsvdcmp(UU,ndatar,ndatar,ndatar,ndatar,w,vv)


! The Picard plot show you eigenvalues. 
if(ran==0) then
   open(66,file=TRIM(output_folder)//'/Picard.out',status='replace')
   do i=1,ndatar
      write(66,*)w(i)
   enddo
   close(66)
endif

!The inversion of Cd was stabilized with truncation of small eigenvalues after singular value decomposition of the system of equations.
!there are no stable analytical formulation for its inverse and determinant. Therefore $\textbf{C}_e^{-1}$ and  $|\textbf{C}_e|$ have to
! be numerically computed with SVD decomposition and removal of a large number of small eigenvalues that destabilize the process.


!Remove small eigenvalues to stabilize the problem

do i=1,ndatar
if (w(i)<wmin) w(i)=0
enddo

!test the svd decomposition ability to solve the equation Cd *x = b
! construct a b 
do i=1,ndatar
 	b(i)=gasdev(ra)
enddo
call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,b,x)
do i=1,ndatar
 	bes(i)=0
         do j=1,ndatar
         	bes(i)=bes(i)+x(j)*Cd(i,j)
         enddo
enddo

like=0
do i=1,ndatar
like=like+real(b(i)-bes(i))**2
enddo
like=like/ndatar

write(*,*)'SVD','test results',sqrt(like),'should be << 1'


!***********************************************************

!                 Get initial likelihood

!***********************************************************


LSr=0
do i=1,ndatar
 br(i)=(d_obsr(i)-wdata(i))
 LSr=LSr+real(br(i)**2)
 enddo
LSr_min=LSr
 br=br/(Ar**2)

! MIsfit INIT with SVD
call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,br,x)

liker=0
do i=1,ndatar
    liker=liker+(d_obsr(i)-wdata(i))*real(x(i))
enddo
liker=liker/(2)
!like=0!
write(*,*)'SVD rank',ran,'like_init',liker

like=liker


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***************************************************************

!    Start the sampling ofthe posterior distribution with the rj_MCMC 

!***************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


av=0
sample=0
th=0
ount=0
PrP=0
PrV=0
PrB=0
PrD=0
AcP=0
AcV=0
AcB=0
AcD=0
post=0
histo=0
histoch=0

ML_Ar=0



best_datar = wdata
like_min = like

do while (sample<nsample)
	ount=ount+1
	
	voro_prop=voro
	like_prop=like
	npt_prop = npt
	Ar_prop = Ar
	AIr_prop = AIr

	u=ran3(ra)
	
	out=1
	move=0
	value=0
	birth=0
	death=0
	noisr=0
	noisd=0
	chwidd=0
	chwidr=0
	logprob=0
	

!*************************************************************

!                   Propose a new model

!*************************************************************


	if (u<0.1) then !change position--------------------------------------------
		move=1
		ind=ceiling(ran3(ra)*npt)
		if (ount.GT.burn_in) then 
			if (voro(ind,1)<(d_max/2)) then
				PrP(1)=PrP(1)+1
			else
				PrP(2)=PrP(2)+1
			endif
		endif

		if (voro(ind,1)<(d_max/2)) then
			 voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd1
		else
			 voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd2
		endif


		!Check if oustide bounds of prior
		if ((voro_prop(ind,1)<d_min).or.(voro_prop(ind,1)>d_max)) then
			out=0
		endif
		call priorvalue(voro_prop(ind,1),d_min,d_max,beta_min,beta_max,width,pv_min,pv_max)
		if ((voro_prop(ind,2)<pv_min).or.(voro_prop(ind,2)>pv_max)) then
			out=0
		endif
	elseif (u<0.3) then ! Change noise parameter for receiver function
		noisr=1
		Prnr = Prnr + 1
		Ar_prop = Ar+gasdev(ra)*pAr
		!Check if oustide bounds of prior
		if ((Ar_prop<Ar_min).or.(Ar_prop>Ar_max)) then
			out=0
		endif

	elseif (u<0.4) then ! change value---------------------------------------------------------
		value=1
		ind=ceiling(ran3(ra)*npt)
		if (ount.GT.burn_in) then 
			if (voro(ind,1)<(d_max/2)) then
				PrV(1)=PrV(1)+1
			else
				PrV(2)=PrV(2)+1
			endif
		endif
		
		if (voro(ind,1)<(d_max/2)) then
			voro_prop(ind,2)=voro(ind,2)+gasdev(ra)*pv1
		else
			voro_prop(ind,2)=voro(ind,2)+gasdev(ra)*pv2
		endif

		!Check if oustide bounds of prior
		call priorvalue(voro_prop(ind,1),d_min,d_max,beta_min,beta_max,width,pv_min,pv_max)
		if ((voro_prop(ind,2)<pv_min).or.(voro_prop(ind,2)>pv_max)) then
			out=0
		endif

	elseif (u<0.7) then !Birth----------------------------------------
			birth = 1
			PrB = PrB + 1
			npt_prop = npt + 1		
			if (npt_prop>npt_max) then
				 out=0
			else
				voro_prop(npt_prop,1) = d_min+ran3(ra)*(d_max-d_min)
				call whichcell(voro_prop(npt_prop,1),voro,npt,npt_max,ind)
				voro_prop(npt_prop,2) = voro(ind,2)+gasdev(ra)*sigmav
				!prob=(1/(sigmav*sqrt(2*pi)))*exp(-(voro(ind,2)-voro_prop(npt_prop,2))**2/(2*sigmav**2))
				logprob=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(npt_prop,2))**2)/(2*sigmav**2)

				!Check bounds					
				call priorvalue(voro_prop(npt_prop,1),d_min,d_max,&
				beta_min,beta_max,width,pv_min,pv_max)
				if ((voro_prop(npt_prop,2)<pv_min).or.(voro_prop(npt_prop,2)>pv_max)) out=0
			
			endif
	else !death!---------------------------------------	
			death = 1
			PrD = PrD + 1
			ind=ceiling(ran3(ra)*npt)
			npt_prop=npt-1
			if (npt_prop<npt_min) then
				 out=0
			else
				voro_prop(ind,:)=voro(npt,:)
				call whichcell(voro(ind,1),voro_prop,npt_prop,npt_max,ind2)
				logprob=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(ind2,2))**2)/(2*sigmav**2)
				!prob=(1/(sigmav*sqrt(2*pi)))*exp(-(voro(ind,2)-voro_prop(ind2,2))**2/(2*sigmav**2))
			endif
	endif!-----------------------------------------------------------------------



!**************************************************************************

!         After  moving a cell, Get the misfit of the proposed model

!**************************************************************************


	if (out==1) then

		! Call Forward MOdel ------------------------------
		call voro2qmodel(voro_prop,npt_prop,npt_max,d_min,d_max,beta,h,vpvs,qa,qb)

!               ooooooooooooooooooooooooooooooooooooooooooooo
		ppara=sin(angle*rad)/v60
		din=asin(ppara*beta(npt_prop)*vpvs(npt_prop))/rad

		call theo(&
		npt_prop, beta, h, vpvs, qa, qb, fs, din,&
		gauss_a, water_c, time_shift, ndatar, wdata )

		!--------compute liker_prop -----------------------
		LSr_prop=0
		do i=1,ndatar
		br(i)=(d_obsr(i)-wdata(i))
		LSr_prop=LSr_prop+br(i)**2
		enddo
		br=br/(Ar_prop**2)

		call svbksb(UU,w,vv,ndatar,ndatar,ndatar,ndatar,br,x)
		
		liker_prop=0
		do i=1,ndatar
		liker_prop=liker_prop+(d_obsr(i)-wdata(i))*x(i)
		enddo
		liker_prop=liker_prop/(2)

		like_prop=liker_prop

		
	endif


!********************************************************************************

!  Now, depending on the type of move, compute acceptance term, and see if we accept or rejetc the model

!********************************************************************************

! See if we accept the proposed change to the model using acceptance ratio.
!   And update the Acceptance ratios, for each type of move. 

	Accept = 0
	if (birth==1)then
        	if (log(ran3(ra))<-log(2*width)-logprob+log(out)-like_prop+like) then
            		accept=1
            		AcB=AcB+1
		endif
        elseif (death==1) then
		if (log(ran3(ra))<log(2*width)+logprob+log(out)-like_prop+like) then
            		accept=1
            		AcD=AcD+1
        	endif

	elseif (noisr==1) then
! Here the likelihood need to be normalised by the determinent of the matrix of datat errors ! lgsigma is the log of the ratio of the determinents of the matrix of data errors.
 		logrsig=ndatar*log(Ar/Ar_prop)
        	if (log(ran3(ra))<logrsig+log(out)-like_prop+like) then
            		accept=1
			Acnr=Acnr+1
        	endif
	
	else !NO JUMP
		if (log(ran3(ra))<log(out)-like_prop+like)then
			accept=1	
             		if((value==1).and.(ount.GT.burn_in)) then
				if (voro(ind,1)<(d_max/2)) then
					AcV(1)=AcV(1)+1
				else
					AcV(2)=AcV(2)+1
				endif
             		elseif((move==1).and.(ount.GT.burn_in))then
				if (voro(ind,1)<(d_max/2)) then
					AcP(1)=AcP(1)+1
				else
					AcP(2)=AcP(2)+1
				endif
			endif
	     	endif!accept
        endif !KIND OF JUMP

	
!***********************************************************************************

!   If we accept the proposed model, update the status of the Markov Chain

!***********************************************************************************

	if (accept==1) then
		voro=voro_prop
		like=like_prop
		liker=liker_prop
		LSr=LSr_prop
		npt=npt_prop
		Ar=Ar_prop
		AIr=AIr_prop	
		if (LSr<LSr_min)then
			 LSr_min = LSr			
			 best_datar=wdata
		
		endif
	endif

!****************************************************************

!                  Store models for ensemble solution

!****************************************************************


IF (ount.GT.burn_in) THEN
	sample=sample+1
	IF (mod(ount,thin)==0) THEN
		th = th + 1
		histo(npt)=histo(npt)+1
		call voro2qmodel(voro,npt,npt_max,d_min,d_max,beta,h,vpvs,qa,qb)
		l=1
		ht=h(l)
		do i=1,disd
			d=(i-1)*prof/real(disd-1)
			if (d<ht)then
				av(i)=av(i)+beta(l)
				mo(i)=beta(l)
			else	
				l=l+1
				av(i)=av(i)+beta(l)
				mo(i)=beta(l)
				if (l<npt) then
					ht=ht+h(l)
				else
					ht=1000
				endif
			endif
			
		enddo
		
		ht=0
		l=1
		ht=h(l)
		do i=1,disd
			d=(i-1)*prof/real(disd-1)
			if (d<ht)then
				v=ceiling((beta(l)-beta_min+width)*&
				disv/(beta_max+2*width-beta_min))
				post(i,v)=post(i,v)+1
		
			else	
				l=l+1
				v=ceiling((beta(l)-beta_min+width)*&
				disv/(beta_max+2*width-beta_min))
				post(i,v)=post(i,v)+1
		
				if (l<npt) then
					ht=ht+h(l)
				else
					ht=1000
				endif
			endif
		enddo
		
		i=ceiling((Ar-Ar_min)*disA/(Ar_max-Ar_min))
		ML_Ar(i) = ML_Ar(i)+1
		
		!Get distribution on changepoint locations.
		ht=0
		do i=1,npt-1
		ht=ht+h(i)
		j=ceiling((ht)*disd/(prof))
		histoch(j)=histoch(j)+1
		enddo	
		
		do i=1,disd
		d=(i-1)*prof/real(disd-1)
		enddo
		
	endif
endif

!get convergence
 conv(ount)=LSr
ncell(ount)=npt
 convAr(ount)=real(Ar)

!**********************************************************************
	
!       Display what is going on every "Display" samples

!**********************************************************************


IF (mod(ount,display).EQ.0) THEN
write(*,*)
write(*,*)'processor number',rank+1,'/',nbproc
write(*,*)'sample:',ount,'/',burn_in+nsample
write(*,*)'number of cells:',npt
if (ount.GT.burn_in) write(*,*)'Acceptance rates'
if (ount.GT.burn_in) write(*,*)'AR_move',100*AcP(1)/PrP(1),100*AcP(2)/PrP(2)
if (ount.GT.burn_in) write(*,*)'AR_value',100*AcV(1)/PrV(1),100*AcV(2)/PrV(2)
write(*,*)'AR_Birth',100*AcB/PrB,'AR_Death',100*AcD/PrD
write(*,*)'AR_Ar',100*Acnr/Prnr
write(*,*)'-----------------------------------------'
write(*,*)
END IF


enddo !End Markov chain


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




av=av/th
!***************************************************************************

! Collect information from all the chains and average everything

!***************************************************************************

call MPI_REDUCE(histoch,histochs,disd,MPI_Integer,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(ML_Ar,ML_Ars,disA,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(av,avs,disd,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(post,posts,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(ncell,ncells,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(conv,convs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(histo,histos,npt_max,MPI_Integer,MPI_Sum,0,MPI_COMM_WORLD,ierror)
call MPI_REDUCE(convAr,convArs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)

avs=avs/nb
 convs=convs/nb
 convArs=convArs/nb
 ncells=ncells/nb



!***********************************************************************

!                      Write the results

!***********************************************************************


IF (ran==0) THEN

open(65,file=TRIM(output_folder)//'/Change_points.out',status='replace')
do i=1,disd
   d=d_min+(i-0.5)*prof/real(disd)
   write(65,*)d,histochs(i)
enddo
close(65)

open(56,file=TRIM(output_folder)//'/Average.out',status='replace')
do i=1,disd
   d=(i-1)*prof/real(disd-1)
   write(56,*)d,avs(i)
enddo
close(56)

open(66,file=TRIM(output_folder)//'/Sigma.out',status='replace')
do i=1,disA
   d=real(Ar_min+(i-0.5)*(Ar_max-Ar_min)/disA)
   write(66,*)d,ML_Ars(i)
enddo
close(66)

open(71,file=TRIM(output_folder)//'/Posterior.out',status='replace')
write(71,*)prof,disd,d_max
write(71,*)beta_min-width,beta_max+width,disv,width
do i=1,disd
   do j=1,disv
      write(71,*)posts(i,j)
   enddo
enddo
close(71)! close the file 


open(72,file=TRIM(output_folder)//'/data_best.out',status='replace')
do i=1,ndatar
   xi=-time_shift+(i-1)/fs
   write(72,*)xi,best_datar(i)
enddo
close(72)


open(54,file=TRIM(output_folder)//'/Convergence_misfit.out',status='replace')
write(54,*)burn_in,nsample
do i=1,nsample+burn_in
   write(54,*)conv(i),convs(i)
enddo
close(54)! close the file 

open(53,file=TRIM(output_folder)//'/Convergence_nb_layers.out',status='replace')
write(53,*)burn_in,nsample
do i=1,nsample+burn_in
   write(53,*)ncell(i),ncells(i)
enddo
close(53)! close the file 

open(52,file=TRIM(output_folder)//'/Convergence_sigma.out',status='replace')
write(52,*)burn_in,nsample
do i=1,nsample+burn_in
   write(52,*)convAr(i),convArs(i)
enddo
close(52)! close the file

open(45,file=TRIM(output_folder)//'/NB_layers.out',status='replace')
do i=1,npt_max
   write(45,*)histos(i)
enddo
close(45)! close the file 

open(46,file=TRIM(output_folder)//'/vpvs.out',status='replace')
do i=1,npt_max
   write(46,*)vpvs(i)
enddo
close(46)! close the file 

endif
CALL cpu_time(t2)
if (ran==0) write(*,*)'time taken by the code was',t2-t1,'seconds'

call MPI_FINALIZE(ierror)! Terminate the parralelization


end! end the main program




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!----------------------------------------------------------------------

!               FUNCTIONS USED BY THE MAIN CODE  

!---------------------------------------------------------------------


!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ----------------------------------------------------------------------------



FUNCTION GASDEV(idum)
implicit none

!     ..Arguments..
integer          idum
real GASDEV
integer iset
!     ..Local..
real v1,v2,r,fac
real ran3

if (idum.lt.0) iset=0
10   v1=2*ran3(idum)-1
v2=2*ran3(idum)-1
r=v1**2+v2**2
if(r.ge.1.or.r.eq.0) GOTO 10
fac=sqrt(-2*log(r)/r)
GASDEV=v2*fac

RETURN
END


!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator
!
! ----------------------------------------------------------------------------

FUNCTION ran3(idum)
implicit none
	
INTEGER idum
INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
! write(*,*)' idum ',idum
if(idum.lt.0.or.iff.eq.0)then
	iff=1
	mj=abs(MSEED-abs(idum))
	mj=mod(mj,MBIG)
	ma(55)=mj
	mk=1
	do 11 i=1,54
	ii=mod(21*i,55)
	ma(ii)=mk
	mk=mj-mk
	if(mk.lt.MZ)mk=mk+MBIG
	mj=ma(ii)
!  write(*,*)' idum av',idum
11      continue
	do 13 k=1,4
	do 12 i=1,55
	ma(i)=ma(i)-ma(1+mod(i+30,55))
	if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
! write(*,*)' idum ap',idum
	inext=0
	inextp=31
	idum=1
endif
! write(*,*)' idum app ',idum
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
!  write(*,*)' idum ',idum
	
return
END

!the code is finished. 
