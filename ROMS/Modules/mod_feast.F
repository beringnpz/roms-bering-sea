#include "cppdefs.h"
      
MODULE mod_feast

#ifdef FEAST

USE mod_param
USE mod_kinds
USE mod_biology

IMPLICIT NONE
     integer, parameter :: MAX_IBM = 15
     integer, parameter :: LENFULL = 1000
     TYPE T_FEAST
			 !id_delM, id_delF, id_delG, id_delRec, id_delProm, id_delCF, id_delCAL
       real(r8), pointer :: ghfal(:,:,:,:,:,:)
       real(r8), pointer :: ghfl(:,:,:,:,:)       
       real(r8), pointer :: ghfsp(:,:,:,:) 
			 real(r8), pointer :: ozm(:,:,:,:)      
			 real(r8), pointer :: zoop_death(:,:,:,:) 
			 real(r8), pointer :: ofdat(:,:,:,:)
			 real(r8), pointer :: incatch(:,:,:,:,:) 
			 !real(r8), pointer :: incatchG(:,:,:)			   
			 real(r8), pointer :: outcatch(:,:,:,:)		  
       real(r8), pointer :: u_swim(:,:,:)
       real(r8), pointer :: v_swim(:,:,:)
       !real(r8), pointer :: eggs(:,:,:)
       real(r8), pointer :: hap_all(:,:,:) 
       real(r8), pointer :: hap_food(:,:,:)
       real(r8), pointer :: hap_pred(:,:,:)  
       real(r8), pointer :: cpy_all(:,:,:)
       real(r8), pointer :: delJ_all(:,:,:)
       real(r8), pointer :: eat_all(:,:)       
       real(r8), pointer :: eat_count(:,:) 
       real(r8), pointer :: mort_all(:,:)
       real(r8), pointer :: eat_simp(:,:)  
			 real(r8), pointer :: eat_plank(:,:)  			  
       real(r8), pointer :: catch_all(:,:,:,:)
       real(r8), pointer :: frate_all(:,:,:,:)
			 
		 ! kya added forcing functions 1/9/2013 
			 real(r8), pointer :: fish_force(:,:,:,:,:)
			 real(r8), pointer :: zoop_force(:,:,:,:,:) 	 
       !delCF(spn) = (W_all(spn) + delW(spn))
       integer, dimension(MAX_IBM) :: sp_ibm, ac_ibm, lc_ibm 
       real(r8), pointer :: PPYY_ibm(:,:,:) 
			 real(r8), pointer :: zoopDC_ibm(:,:,:,:)
       real(r8), pointer :: simpDC_ibm(:,:,:,:)
       real(r8), pointer :: fishDC_ibm(:,:,:,:)
			 real(r8), dimension(LENFULL) :: gCmax,gRfix,gVel_ideal,gRvelV2,gACT,gMzero,gMtwo
			 real(r8), dimension(LENFULL) :: gCpy,gJmax,gCj,gRj,gSFU,gPPYY,gCperV, gEDpy       
     END TYPE T_FEAST
		 
		 TYPE (T_FEAST), allocatable :: GFEAST(:)     
     
     integer, parameter :: id_delM=1,   id_delZ=2,   id_delG=3
		 integer, parameter :: id_delSpawn=4, id_delTZ=5
		 integer, parameter :: id_delCF=6,  id_delCAL=7
		 
		 integer, parameter :: id_forceHappy=1
		 		 
		 		 
		 !real(r8),    parameter :: PREY_ENC_MIN = 1.0e-36
     real(r8),    parameter :: CARBON_TO_WW =  8.83880234228   
     real(r8),    parameter :: WW_TO_CARBON =  0.1131375
!     real(r8),    parameter :: PIPI = 3.14159265358979
     real(r8),    parameter :: G_WW_TO_MG_CARBON = 113.1375
     real(r8),    parameter :: MG_CARBON_TO_G_WW = 0.00883880234228
     integer,     parameter :: NFDAT = 8
 !----------------------------------------------------------------------    
 CONTAINS
 !---------------------------------------------------------------------- 
 SUBROUTINE allocate_gfeast (ng, LBi, UBi, LBj, UBj)

   USE mod_param
   integer, intent(in) :: ng, LBi, UBi, LBj, UBj

   IF (ng.eq.1) allocate ( GFEAST(Ngrids) )

   allocate (GFEAST(ng) % ghfal(LBi:UBi,LBj:UBj,nrates,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES) )
   allocate (GFEAST(ng) % ghfl(LBi:UBi,LBj:UBj,nrates,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS) ) 
   allocate (GFEAST(ng) % ghfsp(LBi:UBi,LBj:UBj,nrates,NUM_SIMPLE_SPECIES) )
   allocate (GFEAST(ng) % ozm(LBi:UBi,LBj:UBj,N(ng),NUM_PLANKTON) )
   allocate (GFEAST(ng) % zoop_death(LBi:UBi,LBj:UBj,N(ng),NUM_PLANKTON) )
	 allocate (GFEAST(ng) % ofdat(LBi:UBi,LBj:UBj,N(ng),NFDAT) )    
	 allocate (GFEAST(ng) % incatch(NUM_GEARS,TOT_FEAST,LBi:UBi,LBj:UBj,2) ) 
	 allocate (GFEAST(ng) % outcatch(LBi:UBi,LBj:UBj,NUM_GEARS,TOT_FEAST) ) 
   allocate (GFEAST(ng) % u_swim(TOT_LINKS,LBi:UBi,LBj:UBj) )
   allocate (GFEAST(ng) % v_swim(TOT_LINKS,LBi:UBi,LBj:UBj) )
   allocate (GFEAST(ng) % hap_all(TOT_LINKS,LBi:UBi,LBj:UBj) )   
   allocate (GFEAST(ng) % hap_pred(TOT_LINKS,LBi:UBi,LBj:UBj) ) 
   allocate (GFEAST(ng) % hap_food(TOT_LINKS,LBi:UBi,LBj:UBj) ) 
   allocate (GFEAST(ng) % cpy_all(TOT_LINKS,LBi:UBi,LBj:UBj) )   
   allocate (GFEAST(ng) % delJ_all(TOT_LINKS,LBi:UBi,LBj:UBj) )  
   !allocate (GFEAST(ng) % eggs(TOT_FEAST,LBi:UBi,LBj:UBj) )
   allocate (GFEAST(ng) % eat_all(TOT_LINKS,TOT_LINKS) )
   allocate (GFEAST(ng) % eat_count(TOT_LINKS,TOT_LINKS) )
   allocate (GFEAST(ng) % mort_all(TOT_LINKS,TOT_LINKS) )
   allocate (GFEAST(ng) % eat_simp(TOT_LINKS,NUM_SIMPLE_SPECIES) )
   allocate (GFEAST(ng) % eat_plank(TOT_LINKS,NUM_PLANKTON) )
	    
   allocate (GFEAST(ng) % fish_force(NFISHFORCE, TOT_FEAST,   LBi:UBi,LBj:UBj,2)) 
   allocate (GFEAST(ng) % zoop_force(NZOOPFORCE, NUM_PLANKTON,LBi:UBi,LBj:UBj,2))    
   
	 allocate (GFEAST(ng) % catch_all(NUM_GEARS,TOT_LINKS,LBi:UBi,LBj:UBj) )
   allocate (GFEAST(ng) % frate_all(NUM_GEARS,TOT_LINKS,LBi:UBi,LBj:UBj) )
   
   allocate (GFEAST(ng) % PPYY_ibm(LBi:UBi,LBj:UBj,MAX_IBM) )
   allocate (GFEAST(ng) % zoopDC_ibm(LBi:UBi,LBj:UBj,MAX_IBM,NUM_PLANKTON) )
   allocate (GFEAST(ng) % simpDC_ibm(LBi:UBi,LBj:UBj,MAX_IBM,NUM_SIMPLE_SPECIES) )
   allocate (GFEAST(ng) % fishDC_ibm(LBi:UBi,LBj:UBj,MAX_IBM,TOT_LINKS) )
   
   RETURN
   END SUBROUTINE allocate_gfeast
 !----------------------------------------------------------------------
   SUBROUTINE initialize_gfeast (ng, tile)

      USE mod_param
      integer, intent(in) :: ng, tile

      integer :: i, j, sp, lc, ac, nr, np
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8

     GFEAST(ng) % sp_ibm = 0
     GFEAST(ng) % ac_ibm = 0
     GFEAST(ng) % lc_ibm = 0
     GFEAST(ng) % PPYY_ibm   = IniVal
		 GFEAST(ng) % zoopDC_ibm = IniVal
     GFEAST(ng) % simpDC_ibm = IniVal
     GFEAST(ng) % fishDC_ibm = IniVal
     GFEAST(ng) % ghfal      = IniVal
     GFEAST(ng) % ghfl       = IniVal
     GFEAST(ng) % ghfsp      = IniVal
     GFEAST(ng) % u_swim     = IniVal
     GFEAST(ng) % v_swim     = IniVal
     GFEAST(ng) % hap_all    = IniVal 
     GFEAST(ng) % cpy_all    = IniVal
     GFEAST(ng) % delJ_all   = IniVal
     !GFEAST(ng) % eggs       = IniVal          
     GFEAST(ng) % incatch    = IniVal
     !GFEAST(ng) % incatch    = IniVal
     GFEAST(ng) % outcatch   = IniVal
     GFEAST(ng) % ozm        = IniVal
     GFEAST(ng) % zoop_death = IniVal
     GFEAST(ng) % ofdat      = IniVal
     GFEAST(ng) % fish_force = IniVal
     GFEAST(ng) % zoop_force = IniVal 
     GFEAST(ng) % eat_all   = IniVal
     GFEAST(ng) % eat_count = IniVal
     GFEAST(ng) % mort_all  = IniVal
		 GFEAST(ng) % eat_simp  = IniVal
		 GFEAST(ng) % eat_plank = IniVal     
     !DO nr = 1,TOT_LINKS
     !   DO np = 1,TOT_LINKS            
     !      GFEAST(ng) % eat_all(nr,np) = IniVal 
     !   end do
     !end do
     

   !print *,tile,"out"  
   END SUBROUTINE initialize_gfeast

!-------------------------------------------------------------------------------
! SUBROUTINE set_mixed(NUM_ROMS_DEPTHS,ROMS_depth,ROMS_temp,ROMS_zoop, &
!                      &mixed, in_T,in_D,zoop,zooptot,ROMS_layer)
!FUNCTION exprate(DailyInc)
!          real(r16), intent(in) :: DailyInc
!					real(r16)             :: exprate
!					real(r16), parameter  :: dt = 144.0
!	exprate = exp(log(MAX(1.0+DailyInc, 1.0e-6))/dt)
!	return
!end function exprate
	
SUBROUTINE set_mixed(LBi,UBi,LBj,UBj,i,j, &
	     & NUM_ROMS_DEPTHS,ROMS_depth,ROMS_edges,ROMS_temp,ROMS_zoop, &
	     & mixed,in_T,in_D,zoop,zooptot,ROMS_layer,ROMS_layer_prop)

     integer, intent(in) :: LBi,UBi,LBj,UBj,i,j 
     integer, intent(in) :: NUM_ROMS_DEPTHS
     real(r8), intent(in) :: ROMS_depth(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS)
     real(r8), intent(in) :: ROMS_edges(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS+1)
	   real(r8), intent(in) :: ROMS_temp(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS)
	   real(r8), intent(in) :: ROMS_zoop(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS,NUM_PLANKTON) 
     real(r8), intent(inout) :: mixed
     real(r8), intent(inout) :: in_T(LAYERS) , in_D(LAYERS)       
	   real(r8), intent(inout) :: zoop(LAYERS,NUM_PLANKTON) !, ZoopMort  ! i,j,d,p  
	   real(r8), intent(inout) :: zooptot(LAYERS,NUM_PLANKTON) !, ZoopMort  ! i,j,d,p  
     integer, intent(inout)  :: ROMS_layer(NUM_ROMS_DEPTHS)  ! i,j,k
     real(r8), intent(inout) :: ROMS_layer_prop(NUM_ROMS_DEPTHS)
     
   integer :: z,l,zz,d,sp,iflag,jflag
   real(r8)     :: MaxT,MinT,MLD,BLD
   integer :: TopLayer, BottomLayer, top, bot
   real(r8)     :: TotVol, DepthVol, cumdep  
	 real(r8)     :: dep(NUM_ROMS_DEPTHS)
	 
	 real(r8), parameter :: FSH_MAX_DEPTH = -300.0

!! HARD LAYER LIMITATION METHOD
   in_T    = 0.0
   in_D    = 0.0
   zoop    = 0.0
   zooptot = 0.0
   
   in_D(1)    = 1.0 	
   in_D(3)    = 1.0
   ROMS_layer = 2
   
    iflag = 99
    jflag = 62 
    
		cumdep = 0.0!cumdep = 0.0
		do d=1,NUM_ROMS_DEPTHS
	    if (ROMS_edges(i,j,d+1) .lt. FSH_MAX_DEPTH) then ! where (g_edges(:,:,d+1) .lt. FEED_DEPTH) 
		     dep(d) = 0.0 !    dep(:,:,d) = 0.0
		  else if (ROMS_edges(i,j,d) .ge. FSH_MAX_DEPTH) then ! elsewhere (g_edges(:,:,d) .ge. FEED_DEPTH) 
		     dep(d) = ROMS_edges(i,j,d+1) - ROMS_edges(i,j,d) !    dep(:,:,d) = g_edges(:,:,d+1) - g_edges(:,:,d)
		  else ! elsewhere
		     dep(d) = ROMS_edges(i,j,d+1) - FSH_MAX_DEPTH !    dep(:,:,d) = g_edges(:,:,d+1) - FEED_DEPTH
		  end if ! end where
		  dep(d) = max(dep(d),0.0)
		  cumdep = cumdep + dep(d) ! cumdep = cumdep + dep(:,:,d)
		  if ((dep(d).gt.0).and.(ROMS_edges(i,j,d+1).gt.ROMS_edges(i,j,d))) then
			   ROMS_layer_prop(d) = dep(d)/(ROMS_edges(i,j,d+1) - ROMS_edges(i,j,d))		  
		  else
		     ROMS_layer_prop(d) = 0.0
			end if
			 
		end do !end do
		do d=1,NUM_ROMS_DEPTHS
		  if (isnan(ROMS_layer_prop(d)) ) then
		    do z=1,NUM_ROMS_DEPTHS
		     write(*,'(a6,3I4,3ES14.4)'),"NDEP",i,j,z,dep(z),ROMS_edges(i,j,z+1),ROMS_edges(i,j,z)
		    end do
		  end if
		end do
		
		in_D(2) = cumdep
		in_T(2) = SUM(ROMS_temp(i,j,:) * dep(:)) / cumdep

		do z = 1,NUM_PLANKTON
		   zooptot(2,z) = MG_CARBON_TO_G_WW * SUM(ROMS_zoop(i,j,:,z) * dep(:))
		end do
		 
!!!!!!!!!!!!!!!!!!!!!!

		!ROMS_temp = ROMS_temp * g_water3d
		!ROMS_avgT = SUM(ROMS_temp(:,:,:) * g_depths(:,:,:), DIM=3) / g_hloc
		!ROMS_depT = SUM(ROMS_temp(:,:,:) *      dep(:,:,:), DIM=3) / cumdep
		!
		!do z= 1,NUM_PLANKTON
    ! ROMS_zoop(:,:,:,z)  = ROMS_zoop(:,:,:,z) * g_water3d
    ! ROMS_sumzoop(:,:,z) = SUM(ROMS_zoop(:,:,:,z) * g_depths(:,:,:), DIM=3) 
		! ROMS_depzoop(:,:,z) = SUM(ROMS_zoop(:,:,:,z) *      dep(:,:,:), DIM=3)    
		!end do

! NO LAYERS				in_T = 0.0    ! in_T(d) = 0.0 !cell_XY->in_T[l]     = 0.0;
! NO LAYERS				in_D = 0.0    ! in_D(d) = 0.0
! NO LAYERS				zoop = 0.0    !zoop(d,sp) = 0.0  !cell_XY->zoop[l][sp] = 0.0;
! NO LAYERS				zooptot = 0.0                
! NO LAYERS   	    MaxT = -100.0 
! NO LAYERS		    MinT =  100.0     
! NO LAYERS     
! NO LAYERS        TotVol = 0.0
! NO LAYERS        DO z = 1,NUM_ROMS_DEPTHS
! NO LAYERS           !  Sum water column volume
! NO LAYERS              TotVol = TotVol + abs(ROMS_depth(i,j,z)) 
! NO LAYERS					 !  Set maximum and minimum temperatures in profile
! NO LAYERS					    if (MaxT<ROMS_temp(i,j,z)) MaxT = ROMS_temp(i,j,z)
! NO LAYERS					    if (MinT>ROMS_temp(i,j,z)) MinT = ROMS_temp(i,j,z)	
! NO LAYERS              
! NO LAYERS        END DO
! NO LAYERS      		
! NO LAYERS   		! First we need to find where the stratification boundaries are
! NO LAYERS   		! using Troy's algorithm
! NO LAYERS        TopLayer    = -1
! NO LAYERS			  BottomLayer = -1    
! NO LAYERS        			
! NO LAYERS        DO zz = NUM_ROMS_DEPTHS, 1, -1
! NO LAYERS           if ((MaxT - ROMS_temp(i,j,zz) > 1.0 ) .AND. (TopLayer<0) ) THEN    
! NO LAYERS			         TopLayer = zz
! NO LAYERS           END IF 
! NO LAYERS        END DO
! NO LAYERS        
! NO LAYERS        DO zz = 1, NUM_ROMS_DEPTHS
! NO LAYERS		       if ((ROMS_temp(i,j,zz) - MinT > 1.0) .AND. (BottomLayer<0)) THEN
! NO LAYERS						  BottomLayer = zz
! NO LAYERS				   END IF      
! NO LAYERS        END DO
! NO LAYERS        
! NO LAYERS
! NO LAYERS      ! Test:  Fixed with stratified 
! NO LAYERS        !TopLayer    = 35
! NO LAYERS        !BottomLayer = 15
! NO LAYERS
! NO LAYERS      ! Test: Fixed with mixed
! NO LAYERS      if (feast_mixed==0) then
! NO LAYERS         TopLayer    = -1
! NO LAYERS         BottomLayer = -1
! NO LAYERS      end if   
! NO LAYERS		  ! Some conditionals for odd water columns (may not work completely)
! NO LAYERS      ! if there's no 1-degree gradient at all (completely mixed)
! NO LAYERS
! NO LAYERS      ! What happens if we force a mixed layer??
! NO LAYERS        !TopLayer    = -1
! NO LAYERS        !BottomLayer = -1
! NO LAYERS      
! NO LAYERS      ! Average and combine depth layers
! NO LAYERS			  if ((TopLayer<0) .OR. (BottomLayer<0)) THEN
! NO LAYERS			      mixed      = 1
! NO LAYERS					  MLD        = 1
! NO LAYERS						BLD        = NUM_ROMS_DEPTHS
! NO LAYERS						in_D(1)  = 1.0 	
! NO LAYERS						in_D(3)  = 1.0 
! NO LAYERS        else
! NO LAYERS					  mixed=0			       
! NO LAYERS					   ! if the layers meet in the middle with no middle
! NO LAYERS					  if (TopLayer<=BottomLayer) THEN				!if (TopLayer>=BottomLayer){
! NO LAYERS							  BottomLayer=TopLayer-1
! NO LAYERS						    if (BottomLayer<1) THEN
! NO LAYERS								  BottomLayer =  1
! NO LAYERS								  TopLayer    =  2
! NO LAYERS							  end if	             !}
! NO LAYERS					   end if			          !}								    
! NO LAYERS					   MLD = TopLayer
! NO LAYERS				 	   BLD = BottomLayer
! NO LAYERS        end if        
! NO LAYERS 
! NO LAYERS      ! Now we loop through layers to set densities and temperatures in each
! NO LAYERS      ! Layer.  If completely mixed, only layer 2 is filled with values.
! NO LAYERS      ! Need to check to make sure odd things don't happen later in this case.
! NO LAYERS 
! NO LAYERS        DO d = 1,LAYERS
! NO LAYERS          if (mixed == 1) then
! NO LAYERS						   top        = NUM_ROMS_DEPTHS
! NO LAYERS						   bot        = 1		
! NO LAYERS					ELSE
! NO LAYERS					 ! Set top and bottom for looping through and Number of Layers
! NO LAYERS						 SELECT CASE(d)
! NO LAYERS								case (1)
! NO LAYERS								      top = NUM_ROMS_DEPTHS 
! NO LAYERS										  bot = TopLayer+1
! NO LAYERS								case (2)  
! NO LAYERS								      top = TopLayer
! NO LAYERS								      bot = BottomLayer+1
! NO LAYERS								case (3)  
! NO LAYERS								      top = BottomLayer
! NO LAYERS											bot = 1
! NO LAYERS								case default
! NO LAYERS							END SELECT
! NO LAYERS         end if					 
! NO LAYERS				 
! NO LAYERS         if (( (mixed == 1) .AND. (d==2) ) .OR. (MIXED == 0)) then 	 
! NO LAYERS			      ! Sum ROMS layers into volume-weighted FEAST layers
! NO LAYERS					    DepthVol=0.0
! NO LAYERS   						DO zz = bot, top 
! NO LAYERS   							 ROMS_layer(zz) = d ! Lookup associating ROMS layer with MLD, etc
! NO LAYERS   							 DepthVol      = DepthVol + ROMS_depth(i,j,zz) 
! NO LAYERS			   				 in_T(d) = in_T(d) + ROMS_temp(i,j,zz) * ROMS_depth(i,j,zz) 
! NO LAYERS								 DO sp = 1,NUM_PLANKTON 
! NO LAYERS                  !zoop(d,sp) = zoop(d,sp) + ROMS_zoop(zz,sp) * ROMS_depth(zz)  
! NO LAYERS                  zooptot(d,sp) = zooptot(d,sp) + ROMS_zoop(i,j,zz,sp) * ROMS_depth(i,j,zz) 
! NO LAYERS                 END DO       
! NO LAYERS              END DO         
! NO LAYERS              if (DepthVol>0.0) then
! NO LAYERS                 in_D(d) = DepthVol
! NO LAYERS                 in_T(d) = in_T(d) / DepthVol
! NO LAYERS							   DO sp = 1,NUM_PLANKTON	
! NO LAYERS                    zoop(d,sp) = zooptot(d,sp) / DepthVol
! NO LAYERS										!zooptot(d,sp) = MG_CARBON_TO_G_WW * zooptot(d,sp);  
! NO LAYERS                 END DO
! NO LAYERS              end if
! NO LAYERS            ! Prevent other division by zero errors later				     
! NO LAYERS				 END IF !}// end else for mixed-layer case
! NO LAYERS		  END DO  !}// end of FEAST layers loopF
! NO LAYERS 					      
! NO LAYERS      DO d = 1,LAYERS
! NO LAYERS				 DO sp = 1,NUM_PLANKTON	 
! NO LAYERS				    zoop(d,sp)    = MG_CARBON_TO_G_WW * zoop(d,sp)    !/ fsh_zoop_wt(sp)
! NO LAYERS				    zooptot(d,sp) = MG_CARBON_TO_G_WW * zooptot(d,sp) !/ fsh_zoop_wt(sp)
! NO LAYERS				 end do
! NO LAYERS			end do
  
 end subroutine set_mixed
 
!--------------------------------------------------------------------------------- 
 
 SUBROUTINE feast_biol(NUM_ROMS_DEPTHS,rmask,ROMS_depth,ROMS_edges,&
        &ROMS_temp,ROMS_zoop,fal,fl,fsp,GF,yday,LBi,UBi,LBj,UBj)

        USE mod_scalars
        USE mod_param
        
        IMPLICIT NONE 
        TYPE (T_FEAST) :: GF
        integer, intent(in) :: LBi,UBi,LBj,UBj,NUM_ROMS_DEPTHS
        real(r8), intent(in) :: ROMS_depth(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS)
        real(r8), intent(in) :: ROMS_edges(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS+1)
	      real(r8), intent(in) :: ROMS_temp(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS)
	      real(r8), intent(in) :: ROMS_zoop(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS,NUM_PLANKTON)
	      real(r8), intent(in) :: rmask(LBi:,LBj:)

        real(r8), intent(inout) :: fal(LBi:UBi,LBj:UBj,nfvaral, &
                               &NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
        real(r8), intent(inout) :: fl(LBi:UBi,LBj:UBj,nfvarl, &
                               &NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS)
        real(r8), intent(inout) :: fsp(LBi:UBi,LBj:UBj,nfvar, &
                               &NUM_SIMPLE_SPECIES)
        real(r8),intent(in) :: yday

      ! Internal indices
        integer  :: i,j,k,z,koff,fvaral,spal,lc,ac,fvarl,spl,fvar,sp
				integer  :: wrap_outflag,ictr,jctr
				
  ! by-layer cell properties
    real(r8) :: mixed
    real(r8), dimension(LAYERS)              :: in_T, in_D      
	  real(r8), dimension(LAYERS,NUM_PLANKTON) :: zoop   
	  real(r8), dimension(LAYERS,NUM_PLANKTON) :: zooptot   	  
    integer, dimension(NUM_ROMS_DEPTHS)      :: ROMS_layer  
    real(r8), dimension(NUM_ROMS_DEPTHS)     :: ROMS_layer_prop 
    
! Run_Step Declarations
  integer :: a,py,lpy,d,ip,pind,age
  real(r8), dimension(LAYERS) :: depth_pref 

              ! Equation 1
                real(r8) TT
				      ! Equation 2 and 2a - weight and biomass
                real(r8) B !,W
              ! Equation 6 - available prey and preference
              ! Preference from equation 7 computed elsewhere and saved
								real(r8) PPYY, zoop_pref !, Smax, Ratio
							! Equations 8 and 10 - diet composition and prey energy density	 
                real(r8) EDpy
                 ! Equation 9 - max. consumption as function of l, temperature
                real(r8) :: Cmax_stom, Cmax_bioe    !double Cmax;
							! Equation 5 // simple for now!
                !real(r8) ACT, Eact 
								!real(r8) :: Ar_W_Br, CC_search, CC_resp, Vel_ideal, Vel_actual								
          			real(r8) ::	pj, ASSIM, CperV, Vc, Fc_T, Cmax, Avel, Jmax
								real(r8) ::	Ar_W_Br, handle, Vr, Fr_T, Rfix, Rvel
								real(r16) :: cubeterm, cA2, cA1, cA0, cQQ, cRR, cDD, cSS, cTT, Vel_ideal
				        real(r8) ::  Cj_assim, Cpy, Cj, RvelV2, Rj, ACT_plus, SFU
							
							! Equations 3, 4, 10	  
                !real(r8) :: TotCons(TOT_FEAST+NUM_PLANKTON)
								real(r8) :: bioen_P, handle_P
                real(r8), dimension(TOT_LINKS) ::  tCpy, tED, tCj, tPP, tHP
                real(r8), dimension(TOT_LINKS,TOT_FEAST+NUM_PLANKTON) :: tTotCons
                ! For forcing vonB growth
                  real(r8) vonB_H, vonB_D, vonB_K, vonB_delW
                ! Equations 11, 12, 13
                  !real(r8) Cj, RFU
                  !real(r8) V, f_T
					      ! SPECIES AND LENGTH Loop II
                  real(r8) recCF												 								               
                  real(r8) Lzero, Lpow, dLen, Ltwo
                  real(r8), dimension(TOT_LINKS) :: delG, delZ, delCAL, delCF, delM, delRec, delSpawn
                  real(r8) propChange, propCAL, propCF, propL
									!real(r8) :: G_len(NUM_AGES), G_lenMinus1(NUM_AGES)
                ! recruitment
                  real(r8) divr, zers
                                              
!  These are values that are needed to pass between loops below.  
! Values needed only within loops are declared within those loops.
   real(r8), dimension(TOT_LINKS) :: mort 
   real(r8), dimension(TOT_LINKS) :: delW, delJ, growJ, spawnJ, totJ 
	 real(r8), dimension(TOT_LINKS) :: fishDC     
   real(r8), dimension(LAYERS,NUM_PLANKTON)    :: ZoopDC
   real(r8), dimension(LAYERS,NUM_PLANKTON)    :: ZoopMort
   real(r8), dimension(LAYERS,NUM_PLANKTON)    :: zoop_WW_eaten
   !real(r8), dimension(LAYERS,NUM_PLANKTON)    :: ZoopMortRate
   real(r8) :: ZoopSurviveRate
   real(r8), dimension(NUM_SIMPLE_SPECIES)     :: SimpDC,SimpleMort,simple_WW_eaten,SimpleMortRate
   real(r8), dimension(NUM_SIMPLE_SPECIES)     :: simple_Bio
   real(r8) :: Ntmp(1:MAX_LCLASS,1:NUM_AGES)
   real(r8) :: N_at_len(MAX_LCLASS)
! Happiness functions
    real(r8) :: happy(LBi:Ubi,LBj:UBj)
    real(r8) :: sad(LBi:Ubi,LBj:UBj)
    real(r8) :: hapXdn(LBi:Ubi,LBj:UBj)
    real(r8) :: hapEdn(LBi:Ubi,LBj:UBj)
    real(r8) :: landXdn(LBi:Ubi,LBj:UBj)
    real(r8) :: landEdn(LBi:Ubi,LBj:UBj)   
    real(r8) :: hapcentx(LBi:Ubi,LBj:UBj)
    real(r8) :: hapcente(LBi:Ubi,LBj:UBj)
    real(r8) :: hapCent(LBi:Ubi,LBj:UBj)
    real(r8) :: hapScaleX(LBi:Ubi,LBj:UBj)
    real(r8) :: hapScaleE(LBi:Ubi,LBj:UBj) 
    real(r8) :: hapx(LBi:Ubi,LBj:UBj) 
    real(r8) :: hape(LBi:Ubi,LBj:UBj) 
   integer :: NN, NS, gr, spn, spy, typeflag, MM,SS,YY,cc, ibm
   integer, dimension(NUM_GEARS*TOT_LINKS) :: gr_sp, gr_gr, gr_link
   integer, dimension(TOT_LINKS) :: murder
   real(r8), dimension(TOT_LINKS) :: N_all,CF_all,W_all,CAL_all  
   real(r8), dimension(NUM_GEARS,TOT_FEAST) :: CPUF, F_required, F_actual
   real real_day
					  
   real_day = yday

wrap_outflag = 0

NN   = ALL_LINKS
GF%ghfal = 0.0
GF%ghfl  = 0.0
GF%ghfsp = 0.0
GF%ozm   = 1.0

GF%eat_all   = 0.0
GF%eat_count = 0.0
GF%mort_all  = 0.0
GF%eat_plank = 0.0
GF%eat_simp  = 0.0
GF%gCmax=0.0; GF%gRfix=0.0; GF%gVel_ideal=0.0; GF%gRvelV2=0.0; GF%gACT=0.0; GF%gMzero=0.0
GF%gPPYY=0.0; GF%gCperV=0.0; GF%gCpy=0.0; GF%gJmax=0.0; GF%gCj=0.0; GF%gRj=0.0; GF%gSFU=0.0
GF%gEDpy=0.0; GF%gMtwo=0.0
  ictr=LBi+(UBi-LBi)/2
  jctr=LBj+(UBj-LBj)/2

 cc=0
! --------------------BEGIN MAIN I,J loop (lasts until end of subroutine)

 
DO i=LBi,UBi
 DO j=LBj,UBj
  
	!if (i .eq. LBi+1 .and. j .eq. LBj+1) then
  !   wrap_outflag=1
  !else
  !   wrap_outflag=0
  !end if
  
  if (rmask(i,j)>0.1) then
  cc=cc+1
! Subroutine to set layers and set temperature and zooplankton concentrations!
! by layers
  call set_mixed(LBi,UBi,LBj,UBj,i,j, &
	     & NUM_ROMS_DEPTHS,ROMS_depth,ROMS_edges,ROMS_temp,ROMS_zoop, &
	     & mixed,in_T,in_D,zoop,zooptot,ROMS_layer,ROMS_layer_prop)
     
	! ZERO THE LOCAL ARRAYS
	     murder = 0				
			 mort = 0.0
			 zoopMort = 0.0
			 zoop_WW_eaten = 0.0
			 simpleMort = 0.0
			 simple_WW_eaten = 0.0
			 N_at_len = 0.0
       !out_catch = 0.0 
       delW      = 0.0
       CPUF = 0.0
       MM = 0       
       simple_bio = 0.0
       
  ! LOAD LOCAL FISH CONDITIONS INTO NN ARRAY
  do sp=1,NUM_SIMPLE_SPECIES
     simple_bio(sp) = fsp(i,j,1,sp)
  end do
  
  do spn=1,NN
   ! sp is the species numbers of aged fish followed by lengthed fish
   ! e.g. 1=POL, 2=COD, 3=ATF, 4=HER, 5=CAP etc.
     sp = sp_all(spn)
   ! Length class (1 through 14 for aged, 1 through 20 for lengthed)
     lc = lc_all(spn)
   ! Age class 1=age 0, 2=age 1, 0 if not an aged fish
     ac = age_all(spn)
          
     !write(4,"(I5) (I5) (I5) (I5)") spn, sp, lc, ac
     if (sp>NUM_AGED_SPECIES) then
        N_all(spn)   = fl(i,j,1,sp - NUM_AGED_SPECIES,lc)
        CF_all(spn)  = fl(i,j,2,sp - NUM_AGED_SPECIES,lc)
        CAL_all(spn) = fl(i,j,3,sp - NUM_AGED_SPECIES,lc)
		 else
        N_all(spn)   = fal(i,j,1,sp,lc,ac)
			  CF_all(spn)  = fal(i,j,2,sp,lc,ac)
        CAL_all(spn) = fal(i,j,3,sp,lc,ac)
				!W_all(spn)  = all_AL_LL_BL(spn) * CF_all(spn)
		 end if
		 if (CF_all(spn)<F_EPSILON) then
		    CF_all(spn) = 1.0
		 end if
     W_all(spn)  = all_AL_LL_BL(spn) * CF_all(spn)
		 totJ(spn)   = W_all(spn) * CAL_all(spn)		 
     do gr = 1, NUM_GEARS
      if (all_QQ(gr,spn)>0.0) then
        if (N_all(spn)> F_EPSILON) then
        !if (fsh_catch_sel(gr,sp)>0.0) then
          !write(4,*) "Instead",NN,i,j,gr,sp,in_catch(i,j,gr,sp)
          if (GF%incatch(gr,sp,i,j,1)>0.0) then
              MM = MM + 1
              gr_sp(MM) = sp
              gr_gr(MM) = gr
              gr_link(MM) = spn
              CPUF(gr,sp)  = CPUF(gr,sp) + all_QQ(gr,spn) * N_all(spn) * W_all(spn)
              if ((i==1) .and. (j==1) .and. (sp==1)) then
							  !write(4,*) "Made it",NN,i,j,gr,sp,lc,ac,MM,CPUF(gr,sp),all_QQ(gr,spn) , N_all(spn), W_all(spn)
							end if
          end if
        end if
      end if
	  end do
		!write(4,*) "All",NN,spn,i,j,sp,lc,ac, labs_all(spn),L_all(spn),N_all(spn), W_all(spn)					  		       
  end do
 	
 	F_required=0.0
 	F_actual=0.0
 	!sum_catch=0.0
  do sp = 1,NUM_AGED_SPECIES+NUM_LENGTHED_SPECIES
	   do gr = 1, NUM_GEARS
        if (CPUF(gr,sp)>0.0) then
	     		 !F_required(gr,sp) = GF%incatch(i,j,gr,sp)/CPUF(gr,sp)
				   !F_actual(gr,sp)   = F_required(gr,sp) !MIN(0.99,F_required(gr,sp))
					 F_actual(gr,sp) = GF%incatch(gr,sp,i,j,1)/CPUF(gr,sp)	 	
        end if
    end do
 end do

  !delF = 0.0
  do spn=1,MM
     gr  = gr_gr(spn)
		 sp  = gr_sp(spn)
		 spy = gr_link(spn)
		 !lc  = lc_all(spy)
		 !ac  = age_all(spy)
		 !write (4,*) 'ha',spn,gr,sp,spy,lc,ac
		 !if (sp>NUM_AGED_SPECIES)  sp = sp - NUM_AGED_SPECIES
		 !delF(spy) = delF(spy) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy)
		 !sum_catch(gr,sp) = sum_catch(gr,sp) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy) * W_all(spy)
     GF%catch_all(gr,spy,i,j) = F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy)
		 GF%frate_all(gr,spy,i,j) = F_actual(gr,sp) * all_QQ(gr,spy)
		     !hfl(i,j,1,sp,lc) = hfl(i,j,1,sp,lc) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy)
		 !    delF(spy) = delF(spy) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy) 
		 !    sum_catch(gr,sp) = sum_catch(gr,sp) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy) * W_all(spy)
		 !else
		 !		 !hfal(i,j,1,sp,lc,ac) = hfal(i,j,1,sp,lc,ac) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy)
		 !		 sum_catch(gr,sp) = sum_catch(gr,sp) + F_actual(gr,sp) * all_QQ(gr,spy) * N_all(spy) * W_all(spy)
		 !		 !write (4,*) 'who',i,j,gr,sp,lc,ac,hfal(i,j,1,sp,lc,ac),F_actual(gr,sp), all_QQ(gr,spy), N_all(spy),W_all(spy)
		 !end if		 
  end do

  !do gr=1,NUM_GEARS
  !   do sp=1,TOT_FEAST         
  !      GFEAST(ng) % outcatch(i,j,gr,sp) = sum_catch(gr,sp)
	!	 end do
	!nd do

  !eggs = 0.0
  !eat_flag = 0
  !if ((i==Istr+1).and.(j==Jstr+1)) then
  !   eat_flag = 1
  !end if
  
  do SS = 1,NN    ! RUN THROUGH IN REVERSE LENGTH ORDER!
     spn = all_order(SS)    
  !do spn = 1,NN
     if (N_all(spn)> -1.) then
		 
		   sp = sp_all(spn)
     
		   ibm=0
       do gr=1,MAX_IBM
        if (sp.eq.GF%sp_ibm(gr) .and. lc_all(spn).eq.GF%lc_ibm(gr) &
			     & .and. age_all(spn).eq.GF%ac_ibm(gr)) ibm=gr
       end do
     ! Zero link-specific quantities
       !PredDC  = 0.0
       fishDC  = 0.0
			 zoopDC  = 0.0
			 simpDC  = 0.0 						        
			 !TotCons = 0.0
     ! Set depth zone to middle if water column structure is completely mixed.
! NO LAYERS				 if (mixed == 1) THEN
					depth_pref(1)=0.0
					depth_pref(2)=1.0
					depth_pref(3)=0.0
! NO LAYERS				 else	
! NO LAYERS						depth_pref(1)=fsh_a_T(sp)
! NO LAYERS						depth_pref(2)=fsh_b_T(sp)
! NO LAYERS						depth_pref(3)=fsh_c_T(sp)
! NO LAYERS				 end if	 
     ! Equation 1
       TT = 0.0
       DO d=1,LAYERS
          TT = TT + depth_pref(d) * in_T(d)
       END DO
              
		! Equation 2 and 2a - weight and biomass
    	B = W_all(spn) * N_all(spn) 
                 
    ! Equation 6 - available prey and preference
    ! Preference from equation 7 computed elsewhere and saved
		! 10/28/11 Change to relative preference
			PPYY = 0.0
		! Non-plankton species
			do YY = SS,NN ! don't have to do this for species longer than this one
				 spy = all_order(YY) !do spy = 1,NN
				 py  = sp_all(spy)
			 ! pre-calculated and stored gamma functions for preferenc
			 !
			 	if (N_all(spy).gt.F_EPSILON) then	
				   fishDC(spy) = N_all(spy) * W_all(spy) * all_preference(spn,spy)
				   PPYY = PPYY + fishDC(spy)
				else
				   fishDC(spy) = 0.0
				end if 
			end do		 					    
	 ! Plankton species
		 DO py = 1,NUM_PLANKTON 
			DO d = 1,LAYERS 
			  if (zooptot(d,py).gt.F_EPSILON) then
				   ZoopDC(d,py) = zooptot(d,py) * depth_pref(d) * zoop_preference(spn,py) 
				   PPYY = PPYY + ZoopDC(d,py)          
				 else
				   ZoopDC(d,py) = 0.0
				 end if                
			END DO
		 END DO								    
	 ! Simple species
		 DO py = 1,NUM_SIMPLE_SPECIES
		   if (simple_Bio(py).gt.F_EPSILON) then
				SimpDC(py) = simple_Bio(py) * simple_preference(spn,py)
				PPYY = PPYY + SimpDC(py)
			 else
			  SimpDC(py) = 0.0
			 end if                         
		 END DO
  
	 ! Equations 8 and 10 - diet composition and prey energy density	 
     if (ibm>0) GF%PPYY_ibm(i,j,ibm) = PPYY
		 EDpy=0.0
     if (PPYY<F_EPSILON) THEN
				EDpy=1.0    
     else            
			  do YY = SS,NN ! don't have to do this for species longer than this one
					 spy         = all_order(YY) !do spy = 1,NN
					 py          = sp_all(spy)
					 fishDC(spy) = fishDC(spy) / PPYY
					 if (ibm>0) GF%fishDC_ibm(i,j,ibm,spy) = fishDC(spy)      
           EDpy        = EDpy + CAL_all(spy) * fishDC(spy)   
				end do   
				DO py = 1,NUM_PLANKTON    
				 DO d = 1,LAYERS   
					 ZoopDC(d,py) = ZoopDC(d,py) / PPYY
					 if (d.eq.2 .and. ibm>0) GF%zoopDC_ibm(i,j,ibm,py) = ZoopDC(d,py)   
					 EDpy         = EDpy + ZoopDC(d,py) * fsh_Zoop_JperG(py)  
				 END DO  
				END DO		
				DO py = 1,NUM_SIMPLE_SPECIES    
					 SimpDC(py) = SimpDC(py) / PPYY
					 if (ibm>0) GF%simpDC_ibm(i,j,ibm,py) = SimpDC(py)  
					 EDpy       = EDpy + SimpDC(py) * fsh_Simple_JperG(py)  
				END DO		
		 end if
		 
		 if (isnan(PPYY)) then
		    do YY = SS,NN 
				 spy = all_order(YY) !do spy = 1,NN
				 py  = sp_all(spy)
				 if (isnan(fishDC(spy))) then
				  write (*,'(A6,6i4,3ES14.4)'),"FISH",i,j,spy,py, &
					& age_all(spy),lc_all(spy), &
					& N_all(spy),W_all(spy),all_preference(spn,spy)			 			  
         end if
        end do
		    DO py = 1,NUM_PLANKTON 
			    DO d = 1,LAYERS 
				   if (isnan(ZoopDC(d,py))) then
					   write (*,'(A6,4i4,3ES14.4)'),"ZOOP",i,j,d,py, &
						 & zooptot(d,py), depth_pref(d), zoop_preference(spn,py)
					 end if                          
			    END DO
		    END DO								    
	 ! Simple species
		    DO py = 1,NUM_SIMPLE_SPECIES 
				  if (isnan(SimpDC(py))) then
				    write (*,'(A6,3i4,3ES14.4)'),"SIMP",i,j,py, &
					  & simple_Bio(py), simple_preference(spn,py) 
					end if                        
		    END DO        
     end if   
  ! Equation 9 - maximum consumption as function of l, temperature
    !Cmax_stom = all_AS_LL_BS(spn) * fsh_A_Q(sp) * exp(fsh_B_Q(sp) * TT)
    !Cmax_bioe = fsh_A_C(sp)*(W_all(spn) ** fsh_B_C(sp) ) * exp(fsh_CQ(sp)*TT)
   
  ! VARIABLE ACTIVITY MODEL Full Excel version
	! Cpy = PPYY * EncRate * (LL_all(spn) ** Enc_Lpow)*(V_act ** Bsearch)
  ! Simplified by assuming Bsearch always equals 1
  ! 13560 Converts g02 to joules  
    !V          = fsh_RTMoverRTM_RT0(sp) - TT / fsh_RTM_RT0(sp)
    !f_T        = (V ** fsh_X_R(sp) ) * exp(fsh_X_R(sp) * (1.0 - V))
    !Ar_W_Br    = fsh_A_R(sp) * (W_all(spn) ** fsh_B_R(sp) ) * 13560.0	
		!CC_search  = PPYY * EDpy * (1.0-fsh_P_j(sp)) * all_Aenc_LL_Benc(spn)
    !CC_resp    = f_T * Ar_W_Br * fsh_A_v(sp) * fsh_B_v(sp)
    !          
		!Vel_ideal  = (CC_search/CC_resp) ** (1.0/(fsh_B_v(sp)-1.0)) 
    !Vel_actual = MAX(Vel_ideal,fsh_minVel(sp))
    !        
    !Cpy        = PPYY * all_Aenc_LL_Benc(spn) * Vel_actual 
		!ACT        = fsh_A_v(sp) * (Vel_actual ** fsh_B_v(sp))  
    !Cj         = Cpy * EDpy
    !RFU        = ACT * f_T * Ar_W_Br  + Cj * fsh_P_j(sp)  
		!				 
		!delJ(spn)           = Cj - RFU
		
	! CUBIC variable activity model (Feb 1 2012)	
        ! LL =StartL
        ! CAL =Cal_inter+Cal_slope*L_all(spn)
        ! CF =1           
        !A_L_B =A_L*L_all(spn)^B_L
        !yearday =1
        !TTT =OFFSET(Sheet2!$D$1,$O5,0)
        !PPYY =OFFSET(Sheet2!$F$1,$O5,0)
        !PREYED =OFFSET(Sheet2!$E$1,$O5,0)
				!PPYY     = PPYY * 100000
				pj       = PPYY * EDpy
				ASSIM    = 1.0 - fsh_P_j(sp)
        !A_W_B = A_L * (L_all(spn) ^ B_L)
        !Wt =A_W_B*CF_all(spn)
        !totJ     = W_all(spn) * CAL_all(spn)
        if (TT.ge.fsh_R_TM(sp)) then
           murder(spn) = 1
           Cpy  = 0.0
			     Cj   = 0.0
			     SFU  = 0.0
           Rj   = 0.0
           Cmax = -2.0
        else if (TT.ge.fsh_C_TM(sp)) then   
           Cpy  = 0.0
			     Cj   = 0.0
			     SFU  = 0.0
			     Cmax = -1.0
           Vr       = fsh_RTMoverRTM_RT0(sp) - TT / fsh_RTM_RT0(sp)
           Fr_T     = (Vr ** fsh_X_R(sp) ) * exp(fsh_X_R(sp) * (1.0 - Vr))
           Ar_W_Br  = fsh_A_R(sp) * (W_all(spn) ** fsh_B_R(sp) ) * 13560.0
					 Rj       = Ar_W_Br * Fr_T
			 else 
        CperV    = PPYY * all_Aenc_LL_Benc(spn) !*EncRate*(L_all(spn)^Enc_Lpow)        
        !Vc       = (CTMj-TT)/(CTMj-CTOj)
         Vc       = fsh_CTMoverCTM_CT0(sp) - TT / fsh_CTM_CT0(sp)
				 Fc_T     = (Vc ** fsh_X_C(sp) ) * exp(fsh_X_C(sp) * (1.0 - Vc))
         Vr       = fsh_RTMoverRTM_RT0(sp) - TT / fsh_RTM_RT0(sp)
         Fr_T     = (Vr ** fsh_X_R(sp) ) * exp(fsh_X_R(sp) * (1.0 - Vr))
         Cmax     = Fc_T * all_AS_LL_BS(spn) !cmax_a * L_all(spn)^cmax_b
         Avel     = EDpy * ASSIM * CperV
         Ar_W_Br  = fsh_A_R(sp) * (W_all(spn) ** fsh_B_R(sp) ) * 13560.0				
				 Rfix     = Ar_W_Br * Fr_T      
         Rvel     = all_Ar_v_LL(spn) * Rfix  !A_v =Av_min+B_Lzero/(1+(B_Lzero/B_Lone-1)*L_all(spn)^B_Lpow)        
				 if (Avel < 1.0) then
				   Vel_ideal = 0.0
				 else
				  Jmax     = EDpy * ASSIM * Cmax
          handle   = 1.0/Jmax
          cubeterm = 2.0*Rvel*handle*handle*Avel*Avel
          cA2      = 4.0*Rvel*handle*Avel/cubeterm
          cA1      = 2.0*Rvel/cubeterm
          cA0      = -Avel/cubeterm
          cQQ       = (3.0*cA1 - cA2*cA2)  / 9.0
          cRR       = (9.0*cA2*cA1 - 27.0*CA0- 2.0*cA2*cA2*cA2) / 54.0
          cDD       = cQQ*cQQ*cQQ + cRR*cRR
          cSS       = (cRR+(SQRT(cDD))) ** (0.3333333333333333333333)!^(1/3)
          cTT       = (cRR-(SQRT(cDD))) ** (0.3333333333333333333333)!^(1/3)
          Vel_ideal = -(cA2/3.0) + (cSS+cTT)
         end if				  
				Cj_assim  = EDpy * ASSIM * CperV * Vel_ideal/ &
				          & (1.0 + handle * EDpy * ASSIM * CperV * Vel_ideal)
        Cpy       = Cj_assim / (EDpy*ASSIM)
        Cj        = Cpy * EDpy
        RvelV2    = Rvel*Vel_ideal*Vel_ideal
        Rj        = Rfix + RvelV2
        ACT_plus  = RvelV2/Rfix
        SFU       = fsh_P_j(sp) * Cj
			end if
			 
      ! recruitment
        delSpawn(spn) = 0.0
        if (yday.ge.fsh_sp_sday(sp) .and. yday.le.fsh_sp_eday(sp) ) then
             divr = (fsh_sp_eday(sp) - fsh_sp_sday(sp) + 1)
             delSpawn(spn)   = W_all(spn) * all_mat_rec(spn) / divr
        end if
        
      ! Summing the joules  
        growJ(spn)  = Cj - Rj - SFU
        spawnJ(spn) = delSpawn(spn) * CAL_all(spn)
      
      if ((i.eq.ictr).and.(j.eq.jctr)) then
        GF%gPPYY(spn)      = PPYY
        GF%gEDpy(spn)      = EDpy
        GF%gCmax(spn)      = Cmax
        GF%gCperV(spn)     = CperV
				GF%gRfix(spn)      = Rfix
				GF%gVel_ideal(spn) = Vel_ideal
				GF%gRvelV2(spn)    = RvelV2
				GF%gACT(spn)       = 1.0 + ACT_plus
				GF%gCpy(spn)       = Cpy
				GF%gJmax(spn)      = Jmax
				GF%gCj(spn)        = Cj
				GF%gRj(spn)        = Rj
				GF%gSFU(spn)       = SFU
      end if
      
        !if (wrap_outflag .eq. 1 .and. sp .eq. 1) then
        if (isnan(Cpy)) then
				   write(*,'(A10,i6)'),"FDNAN",spn
           write(*,'(A10,I6,G14.6)'),"PPYY",spn,PPYY
           write(*,'(A10,I6,G14.6)'),"EDpy",spn,EDpy
           write(*,'(A10,I6,G14.6)'),"TT",spn,TT
           write(*,'(A10,I6,G14.6)'),"LL",spn,L_all(spn)
           write(*,'(A10,I6,G14.6)'),"CF",spn,CF_all(spn)
           write(*,'(A10,I6,G14.6)'),"CAL",spn,CAL_all(spn)           
           write(*,'(A10,I6,G14.6)'),"Wt",spn,W_all(spn)
           write(*,'(A10,I6,G14.6)'),"Aenc_LL_B",spn,all_Aenc_LL_Benc(spn) 
           write(*,'(A10,I6,G14.6)'),"CperV",spn,CperV           
           write(*,'(A10,I6,G14.6)'),"Avel",spn,Avel           
           write(*,'(A10,I6,G14.6)'),"Cmax",spn,Cmax           
           write(*,'(A10,I6,G14.6)'),"Jmax",spn,Jmax 
           write(*,'(A10,I6,G14.6)'),"handle",spn,handle 
           write(*,'(A10,I6,G14.6)'),"Fr_T",spn,Fr_T          
           write(*,'(A10,I6,G14.6)'),"Rfix",spn,Rfix
           write(*,'(A10,I6,G14.6)'),"Av",spn,all_Ar_v_LL(spn)				   
           write(*,'(A10,I6,G14.6)'),"Rvel",spn,Rvel
					 write(*,'(A10,I6,G14.6)'),"cubeterm",spn,cubeterm
					 write(*,'(A10,I6,G14.6)'),"cA2",spn,cA2
					   write(*,'(A10,I6,G14.6)'),"cA1",spn,cA1
					   write(*,'(A10,I6,G14.6)'),"cA0",spn,cA0
					   write(*,'(A10,I6,G14.6)'),"cQQ",spn,cQQ
					   write(*,'(A10,I6,G14.6)'),"cRR",spn,cRR
					   write(*,'(A10,I6,G14.6)'),"cDD",spn,cDD
					   write(*,'(A10,I6,G14.6)'),"cSS",spn,cSS
					   write(*,'(A10,I6,G14.6)'),"cTT",spn,cTT					   
					 write(*,'(A10,I6,G14.6)'),"Vel_ideal",spn,Vel_ideal
					 write(*,'(A10,I6,G14.6)'),"Cpy",spn,Cpy					 				 				 					            
        end if
        !NewLL     = ((NewJ/CAL_all(spn))/A_L)^(1/B_L)
        !CJ/J =Cj/totJ
        !RJ/j =(Rj+SFU)/totJ
        !gro_prop  = growJ/totJ
        !FitProp =(I6-I5)/I5
        !Cmax =Cmax/W_all(spn)
        !Cjassim/Cmaxassim =Cj_assim/Jmax
						  
		! 10/6/11 temporary delW for purposes of happiness
		  !delW(spn) = (Cj - RFU)/all_ED(spn)

  ! No more changes to Cpy, it should be set at this point
    !bioen_P = Cpy / Cmax_bioe

! NOOUT    if ((i .eq. (88-1)).and.(j .eq. (60-1))) then
! NOOUT		 if (sp .eq. 1) then
! NOOUT			if (lc_all(spn) .eq. 7) then
! NOOUT			   if (age_all(spn) .eq. 1) then
! NOOUT			    do z = 1,NUM_PLANKTON
! NOOUT			     do k = 1,NUM_ROMS_DEPTHS
! NOOUT			        write(*,'(4I4,2G18.8)'),i,j,k,z, &
! NOOUT							& ROMS_zoop(i,j,k,z), ROMS_depth(i,j,k)
! NOOUT						end do	
! NOOUT			    end do                                               
! NOOUT			   end if
! NOOUT				 write(*,'(4I6,F6.2)'),i,j,age_all(spn),lc_all(spn),real(L_all(spn))
! NOOUT			   write(*,'(A16,G18.8)'),"N_all",real(N_all(spn))
! NOOUT			   write(*,'(A16,G18.8)'),"CF_all",real(CF_all(spn))
! NOOUT			   write(*,'(A16,G18.8)'),"W_all",real(W_all(spn))
! NOOUT         write(*,'(A16,G18.8)'),"all_ED",real(all_ED(spn))
! NOOUT         write(*,'(A16,G18.8)'),"PPYY",real(PPYY)
! NOOUT         write(*,'(A16,G18.8)'),"EDpy",real(EDpy)
! NOOUT         write(*,'(A16,G18.8)'),"TT",real(TT)
! NOOUT         !write(*,'(A16,G14.6)'),"V",real(V)
! NOOUT         !write(*,'(A16,G14.6)'),"X_R",real(fsh_X_R(sp))
! NOOUT      	 !write(*,'(A16,G14.6)'),"RTM-RT0",real(fsh_RTM_RT0(sp))
! NOOUT	  		 !write(*,'(A16,G14.6)'),"over",real(fsh_RTMoverRTM_RT0(sp))
! NOOUT				 write(*,'(A16,G18.8)'),"f_T",real(f_T)
! NOOUT         write(*,'(A16,G18.8)'),"Ar_W_Br", real(Ar_W_Br)
! NOOUT         write(*,'(A16,G18.8)'),"CC_search",real(CC_search)
! NOOUT         write(*,'(A16,G18.8)'),"CC_resp",real(CC_resp)
! NOOUT         write(*,'(A16,G18.8)'),"Vel_ideal",real(Vel_ideal)
! NOOUT         write(*,'(A16,G18.8)'),"Vel_actual",real(Vel_actual)
! NOOUT         write(*,'(A16,G18.8)'),"Cpy",real(Cpy)
! NOOUT         write(*,'(A16,G18.8)'),"ACT",real(ACT)
! NOOUT         write(*,'(A16,G18.8)'),"Cj",real(Cj)
! NOOUT         write(*,'(A16,G18.8)'),"RFU",real(RFU)
! NOOUT         write(*,'(A16,G18.8)'),"delJ",real(delJ(spn))
! NOOUT         write(*,'(A4,A6,A4,A4,4A14)')," PdL","Spec","A_D","PyL",&
! NOOUT				      & "BioPy_gWWm2","Pref","DC","PPYYtot"       
! NOOUT				 DO d = 1,LAYERS   
! NOOUT					DO py = 1,NUM_PLANKTON 
! NOOUT						 if (ZoopDC(d,py) .gt. 1e-20) then
! NOOUT							write(*,'(F4.0,A4,I2.2,I4,F4.0,4G14.6)'),real(L_all(spn)), &
! NOOUT							& 'ZZ_',py,d,0.0, &
! NOOUT							& real(zooptot(d,py) * depth_pref(d)), &
! NOOUT							& real(zoop_preference(spn,py)),real(ZoopDC(d,py)), &
! NOOUT							& real(zooptot(d,py) * depth_pref(d) * zoop_preference(spn,py))
! NOOUT						 end if 
! NOOUT					END DO
! NOOUT				 END DO
! NOOUT				 DO py = 1,NUM_SIMPLE_SPECIES    
! NOOUT						 if (simpDC(py) .gt. 1e-20) then
! NOOUT							write(*,'(F4.0,A4,I2.2,I4,F4.0,4G14.6)'),real(L_all(spn)), &
! NOOUT							& 'SS_',py, 0, 0.0, &
! NOOUT							& real(simple_Bio(py)), &
! NOOUT							& real(simple_preference(spn,py)),real(simpDC(py)), &
! NOOUT							& real(simple_Bio(py) * simple_preference(spn,py))												                           
! NOOUT						 end if
! NOOUT				END DO
! NOOUT				do YY = SS,NN
! NOOUT						spy = all_order(YY) !do spy = 1,NN
! NOOUT						 if (fishDC(spy) .gt. 1e-20) then
! NOOUT							write(*,'(F4.0,A4,I2.2,I4,F4.0,4G14.6)'),real(L_all(spn)), &
! NOOUT							&	'FF_',sp_all(spy), age_all(spy), L_all(spy), &
! NOOUT							& real(N_all(spy) * W_all(spy)), &
! NOOUT							& real(all_preference(spn,spy)),real(fishDC(spy)), &
! NOOUT							& real(N_all(spy) * W_all(spy) * all_preference(spn,spy) )												                           
! NOOUT						 end if
! NOOUT			 END DO
! NOOUT			end if
! NOOUT		 end if
! NOOUT    end if
								
								
       ! Equation 14 predation mortality
				 do YY = SS,NN ! don't have to do this for species longer than this one
						spy=all_order(YY) !do spy = 1,NN                             
            mort(spy) = mort(spy) + N_all(spn) * Cpy * fishDC(spy)/W_all(spy)
            !TotCons(sp_all(spy)) = TotCons(sp_all(spy))+N_all(spn)*Cpy*fishDC(spy)
            if ((i.eq.ictr).and.(j.eq.jctr)) then
             GF % eat_all(spn,spy)    = Cpy * fishDC(spy)
					   GF % eat_count(spn,spy)  = Cpy * fishDC(spy)/W_all(spy)
					   if (N_all(spy).gt.0) then
						   GF % mort_all(spy,spn)   =  365.0 * log(max(1e-6, 1.0 - &
						   & ((N_all(spn) * Cpy * fishDC(spy)/W_all(spy))/(N_all(spy)))))  
						 else
						   GF % mort_all(spy,spn)   =  0.0
						 end if
						end if
         END DO    
								             
        ! This conversion is for trawl survey data per hectare
          ! Units:  DC * g/pred * pred/10,000m^2(1 hectare) * 1/depth(m) * 1000mg/g  = mg/m^3
        ! This one is for pred/m^2
          ! Units:  DC * g/pred * pred/m^2(1 hectare) * 1/depth(m) * 1000mg/g  = mg/m^3
          DO py = 1,NUM_PLANKTON    !for (py=0; py<NUM_PLANKTON; py++){
             DO d = 1, LAYERS
                 if (in_D(d)>0) THEN
                     zoop_WW_eaten(d,py) = zoop_WW_eaten(d,py) + ZoopDC(d,py) * Cpy * N_all(spn) 
                     if ((i.eq.ictr).and.(j.eq.jctr)) then
                         GF%eat_plank(spn,py) = GF%eat_plank(spn,py) + ZoopDC(d,py)*Cpy
										 end if										 
								 else
               	     zoop_WW_eaten(d,py) = 0.0
								 end if												
		   					 !TotCons(TOT_FEAST+py) = TotCons(TOT_FEAST+py) + ZoopDC(d,py) * Cpy * N_all(spn)

             END DO								
         END DO				

        ! This loop is for the simple (not-lenghted) FEAST species
          DO py = 1,NUM_SIMPLE_SPECIES    !for (py=0; py<NUM_PLANKTON; py++){
             simple_WW_eaten(py) = simple_WW_eaten(py) + SimpDC(py) * Cpy * N_all(spn)
             if ((i.eq.ictr).and.(j.eq.jctr)) then
                 GF%eat_simp(spn,py) = GF%eat_simp(spn,py) + SimpDC(py)*Cpy
						 end if							  												
	 	    		 !TotCons(TOT_LENGTHED+py) = TotCons(TOT_LENGTHED+py) + SimpDC(py) * Cpy * N_all(spn)					
         END DO	
     	    
							 tCpy(spn) = Cpy
							 tED(spn)  = EDpy
							 tCj(spn)  = Cj
							 tPP(spn)  = bioen_P
							 tHP(spn)  = handle_P
							 
							 !do py = 1,TOT_FEAST+NUM_PLANKTON
							 !   tTotCons(spn,py) = TotCons(py)
							 !end do 


        end if  !end of if N>0 loop
      
	END DO!} //end of first species spn loop -------------------------------------

  delCF = 0.0
  delG  = 0.0
  delM  = 0.0
  !delF  = 0.0
  
  do SS = 1,NN    ! RUN THROUGH IN REVERSE LENGTH ORDER!
     spn = all_order(SS)
		 sp  = sp_all(spn)
     lc = lc_all(spn)
     ac = age_all(spn)   
     !gg = 0.0
            ! Equations 17, 18
           ! Not allowed to get shorter, just thinner
             !if (delW(spn)>0.0)     gg=fsh_g_W(sp) ! else 0
	    			 !if (next_len(spn)==0)  gg=0.0
             !delL = (gg*delW(spn)/(CF_all(spn)*fsh_A_L(sp) ) + (L_all(spn)**fsh_B_L(sp) )) ** (1.0/fsh_B_L(sp) ) 
             !delL  = delL - L_all(spn)
             !delCF(spn) = (W_all(spn) + delW(spn))/((fsh_A_L(sp) *((delL+L_all(spn))**fsh_B_L(sp) )))
             !delCF(spn) = delCF(spn) - CF_all(spn)             
     
        !delJ(spn) = newJ - totJ				
				!delW(spn) = delW(spn) - delSpawn(spn)
												
			! Translate joules to state variables	
        ! Prop_dense =IF(PropChange<1,PropChange^(2/3),IF(CAL_all(spn)<MaxDen+Cal_0+L_all(spn)*(Cal_80-Cal_0)/80,MIN(1.01,SQRT(PropChange)),1))
        ! Prop_CF =IF(PropChange<1,PropChange^(1/3),IF(CF_all(spn)<MaxCF,MIN(1.01,SQRT(PropChange)),1))
        ! Prop_W =PropChange/(Prop_dense*Prop_CF)

     propChange = MAX(totJ(spn)+growJ(spn)-spawnJ(spn),0.01) / totJ(spn)
		 propCAL    = 1.0
		 propCF     = 1.0
		 if ((propChange .le. 1.0).or.(next_len(spn).eq.0)) then
				propCAL = propChange ** (0.6666666666666666666667)
				propCF  = propChange ** (0.3333333333333333333333)
				propL   = 1.0
				dLen    = 0.0
				delG(spn)  = 0.0
		 else
		    ! 7/18/2012 Removed Minimum 1.01 change
				!if (CF_all(spn) < 1.05) propCF = MIN(1.01,SQRT(PropChange))
				!if (CAL_all(spn)< 1.05* CALbase_all(spn)) propCAL = MIN(1.01,SQRT(PropChange))
				if (CF_all(spn) < 1.05) propCF = SQRT(PropChange)
				if (CAL_all(spn)< 1.05* CALbase_all(spn)) propCAL = SQRT(PropChange)
		    propL      = propChange/(propCAL*propCF)
        dLen       = (propL * all_AL_LL_BL(spn)/fsh_A_L(sp))**(1/fsh_B_L(sp)) - L_all(spn)
        delG(spn)  = MAX(0.0,MIN(0.99,MAX(0.0,dLen/L_width(spn))))		    		    
		 end if
		                   
		! Change Accounting	NOW AS RATES	       
      delCAL(spn) = propCAL - 1.0 ! * CAL_all(spn)
      delCF(spn)  = propCF  - 1.0 ! * CF_all(spn)
				  
      !if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
      !   if ((sp .eq. 1)) then     
			!			write (*,'(A6,F14.2,4I4,ES14.4)'),"FSPW",yday,i,j,ac,lc, &
			!	           & delSpawn(spn)          
			!			!write (*,'(A6,F14.2,4I4,9ES14.4)'),"FGRO",yday,i,j,ac,lc, &
			!	    !       & totJ(spn),growJ(spn),propCAL,propCF,propL, &
			!			!			 & delCAL(spn),delCF(spn),dLen,delG(spn)
			!	 end if
			!end if
     
    ! Mortality, equation 15, 14
      Lpow = (fsh_mu(sp) * CF_all(spn))** fsh_zeta(sp)
      Lzero = fsh_omega(sp)/365.0 * (Lpow + fsh_smallomega(sp) )/Lpow
      !if ((i.eq.ictr).and.(j.eq.jctr)) then
			!   GF%gMzero(spn) = -365.0 * Lzero  
			!end if
    ! convert yearly mortality rate to daily loss, might
    ! make function of survival not mort to save an exp() ?			
      
		! KYA 5/29/12 was this in error?  Down there it's a rate	
			! Lzero = (1.0 - exp(-Lzero)) 

    ! Natural Mortality (INCLUDES PREDATION AND Mzero)
      if (N_all(spn) .gt. 0.0) then      
				 !delM(spn) = MAX(0.0,mort(spn)/N_all(spn))
				 Ltwo      = MAX(1e-6,(N_all(spn)-mort(spn))/N_all(spn))
				 delM(spn) = -log(Ltwo) !MAX(0.0,mort(spn)/N_all(spn))
				 delZ(spn) = MAX(0.0,Lzero)
      else
         delM(spn) = 0.0
         delZ(spn) = 0.0
      end if

      if ((i.eq.ictr).and.(j.eq.jctr)) then
			   GF%gMzero(spn) = -365.0 * delZ(spn)  
			   GF%gMtwo(spn)  = -365.0 * delM(spn) 
			end if
                                       
   END DO !}   // end of second lengthclass loop -------------------------------
    
  ! OUTPUT RATE VARIABLES FOR TIME STEPPING
    do SS = 1,NN    ! RUN THROUGH IN REVERSE LENGTH ORDER!
       spn = all_order(SS)
		   sp  = sp_all(spn)
       lc = lc_all(spn)
       ac = age_all(spn) 
         ! KYA 10/9/11 replace delN here with individual rates
				 ! delN(spn) = -delM(spn)-delF(spn)-delG(spn)+delRec(spn)
         ! if (prev_len(spn)>1)  delN(spn) = delN(spn) + delG(prev_len(spn))         
       
       ! KYA 5/29/12 HAPPINESS OPTIONS
			 ! Growth only
			 !  GF%hap_all(spn,i,j)  = growJ(spn)/totJ(spn)
       ! All Biomass
         !GF%hap_all(spn,i,j)  = growJ(spn)/totJ(spn) + (1.0-N_all(spn)*exp(-(delM(spn) + delZ(spn))))
         !BJchange_fix <- ((LOWLIM+hdat$PdNN) * hdat$GrowJ  +  hdat$TotJ * Nchange)/((LOWLIM+hdat$PdNN) * hdat$TotJ)

              
      ! KYA 03-11-15 COMMENTED OUT FOLLOWING LINE to split happiness into pieces 
      !  GF%hap_all(spn,i,j) = ((F_EPSILON+N_all(spn))*growJ(spn)  +     &
			!	   &  totJ(spn)*N_all(spn)*(exp(-delM(spn)-delZ(spn)) - 1.0)) / &
			!		 &  ((F_EPSILON+N_all(spn))*totJ(spn))
			
      ! Units at this point should be % net population growth (in calories) 
      ! and % of population calories surviving predation	
        GF%hap_food(spn,i,j) = growJ(spn) * (F_EPSILON+N_all(spn)) / &
                             & ((F_EPSILON+N_all(spn))*totJ(spn))  
        GF%hap_pred(spn,i,j) = (totJ(spn) * N_all(spn) * &
                             & (exp(-delM(spn)-delZ(spn)) - 1.0))/ & 
                               ((F_EPSILON+N_all(spn))*totJ(spn))  	        
        
        GF%hap_all(spn,i,j)  = GF%hap_food(spn,i,j) + GF%hap_pred(spn,i,j)	 
        
        
        	 
		   GF%cpy_all(spn,i,j)  = tCpy(spn)/W_all(spn)
		   GF%delJ_all(spn,i,j) = delJ(spn)   
 			 			 
			 !GF%eggs(sp,i,j)      = GF%eggs(sp,i,j) +&
       !                     & N_all(spn) * delSpawn(spn) * fsh_rec_prop(sp)
       if (CF_all(spn).lt.(0.25))                   murder(spn) = 1
			 if (CAL_all(spn).lt.(0.25*CALbase_all(spn))) murder(spn) = 1                     
                            
       if (sp > NUM_AGED_SPECIES) then
          sp = sp - NUM_AGED_SPECIES ! NOTE NOAGE LENGTHED SPECIES ARE FIRST
          GF%ghfl(i,j,id_delSpawn,sp,lc)     = delSpawn(spn) * fsh_rec_prop(sp)/144.0
					GF%ghfl(i,j,id_delM,sp,lc)         = exp(-delM(spn)/144.0)             ! RATE
          GF%ghfl(i,j,id_delZ,sp,lc)         = exp(-delZ(spn)/144.0)             ! RATE
          GF%ghfl(i,j,id_delTZ,sp,lc)        = exp((-delZ(spn)-delM(spn))/144.0)
					if (murder(spn).gt.0) GF%ghfl(i,j,id_delTZ,sp,lc) = 0.0
					GF%ghfl(i,j,id_delG,sp,lc)         = exp(log(max(1e-6,1.0 -   delG(spn)))/144.0)   ! RATE
				 	GF%ghfl(i,j,id_delCF,sp,lc)        = exp(log(max(1e-6,1.0 +  delCF(spn)))/144.0)  ! RATE
          GF%ghfl(i,j,id_delCAL,sp,lc)       = exp(log(max(1e-6,1.0 + delCAL(spn)))/144.0)  ! RATE				
       else
          GF%ghfal(i,j,id_delSpawn,sp,lc,ac) = delSpawn(spn) * fsh_rec_prop(sp)/144.0
          GF%ghfal(i,j,id_delM,sp,lc,ac)     = exp(-delM(spn)/144.0)    ! RATE
          GF%ghfal(i,j,id_delZ,sp,lc,ac)     = exp(-delZ(spn)/144.0)   ! RATE
          GF%ghfal(i,j,id_delTZ,sp,lc,ac)    = exp((-delZ(spn)-delM(spn))/144.0) ! RATE
	        if (murder(spn).gt.0) GF%ghfal(i,j,id_delTZ,sp,lc,ac) = 0.0
					GF%ghfal(i,j,id_delG,sp,lc,ac)     = exp(log(max(1e-6,1.0 -   delG(spn)))/144.0)   ! RATE
					GF%ghfal(i,j,id_delCF,sp,lc,ac)    = exp(log(max(1e-6,1.0 +  delCF(spn)))/144.0)  ! RATE
          GF%ghfal(i,j,id_delCAL,sp,lc,ac)   = exp(log(max(1e-6,1.0 + delCAL(spn)))/144.0)  ! RATE
					
      !if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
      !   if ((sp .eq. 1)) then     
      !      write (*,'(A6,F14.2,4I4,7ES18.8)'),"FRATE",yday,i,j,ac,lc, &
		  !		           & GF%ghfal(i,j,id_delSpawn,sp,lc,ac), &
			!						 & GF%ghfal(i,j,id_delM,sp,lc,ac),  &
			!	           & GF%ghfal(i,j,id_delZ,sp,lc,ac),  &
			!						 & GF%ghfal(i,j,id_delTZ,sp,lc,ac), &
			!						 & GF%ghfal(i,j,id_delG,sp,lc,ac),  &
			!						 & GF%ghfal(i,j,id_delCF,sp,lc,ac), &
			!						 & GF%ghfal(i,j,id_delCAL,sp,lc,ac)
			!	 end if
			!end if					
					     
       end if
 
    end do	

! VERTICAL VERSION		    DO py = 1,NUM_PLANKTON
! VERTICAL VERSION		       DO d = 1,LAYERS							 
! VERTICAL VERSION             	if (zooptot(d,py)>0 .AND. in_D(d)>0) THEN						
! VERTICAL VERSION							    ZoopMortRate(d,py) = zoop_WW_eaten(d,py) / (zooptot(d,py))
! VERTICAL VERSION							else
! VERTICAL VERSION							    ZoopMortRate(d,py) = 0.0
! VERTICAL VERSION							end if
! VERTICAL VERSION							
! VERTICAL VERSION           END DO ! End of Layers loop (top, middle, bottom ocean layers)
! VERTICAL VERSION				END DO	! end of Plankton py loop
! VERTICAL VERSION
! VERTICAL VERSION        DO z=1, NUM_ROMS_DEPTHS
! VERTICAL VERSION           DO sp=1,NUM_PLANKTON
! VERTICAL VERSION             !if (feast_coupled>0) then
! VERTICAL VERSION                GF % ozm(i,j,z,sp) = ROMS_layer_prop(z) * &
! VERTICAL VERSION								   & ROMS_zoop(i,j,z,sp) * ZoopMortRate(ROMS_layer(z),sp)
! VERTICAL VERSION		       END DO
! VERTICAL VERSION        END DO

! All in Layer 2 version
  d = 2
  DO py = 1,NUM_PLANKTON
   if (zooptot(d,py)>0 .AND. in_D(d)>0) THEN						
    DO z=1, NUM_ROMS_DEPTHS
     ZoopSurviveRate  = min(1.0, max(1e-6, 1.0-ROMS_layer_prop(z)*zoop_WW_eaten(d,py)/(zooptot(d,py))))	  
     GF%ozm(i,j,z,py) = exp(log(ZoopSurviveRate)/144.0)
	   if ((i.eq.(121)).and.(j.eq.51)) then
	      write(*,'(A6,2I4,6ES14.4)'),"ZMORT",py,z,zooptot(d,py),in_D(d), &
				& ROMS_layer_prop(z),zoop_WW_eaten(d,py),ZoopSurviveRate, &
				& GF%ozm(i,j,z,py)
		 end if
		END DO
	 else
     GF % ozm(i,j,:,py) = 1.0
   end if
  END DO   

    end if ! END OF RMASK IF
  END DO ! main j loop
END DO  ! main i loop
! //*********************************************************************//

! 11/6/11 HAPPINESS MOVED HERE FROM FEAST_STEP
do spn=1,ALL_LINKS
   sp = sp_all(spn)
   lc = lc_all(spn)
   ac = age_all(spn)
          
   happy = 1./(1.+exp(-swim_selk(spn)*(GF%hap_all(spn,:,:)-swim_selm(spn))))
   happy = happy * (1. - GF%fish_force(1,sp_all(spn),:,:,1))
   sad    = 1.0 - happy
    hapXdn  = 0.0
    hapEdn  = 0.0
    landXdn = 0.0
    landEdn = 0.0
               
    hapXdn((LBi+1):UBi,:) = happy(LBi:(UBi-1),:)
    hapEdn(:,(LBj+1):UBj) = happy(:,LBj:(UBj-1)) 
   
    landXdn((LBi+1):UBi,:) = rmask(LBi:(UBi-1),:)
    landEdn(:,(LBj+1):UBj) = rmask(:,LBj:(UBj-1))   

    hapx   = (happy - hapXdn) * rmask * landXdn
    hape   = (happy - hapEdn) * rmask * landEdn
        
    hapcentx = 0.0
    hapcente = 0.0
    do i = (LBi+1),(UBi-1)
     do j = (LBj+1),(UBj-1)
      hapcentx(i,j)=0.25*(hapx(i,j)+hapx(i+1,j)+hapx(i,j-1)+hapx(i+1,j-1))
      hapcente(i,j)=0.25*(hape(i,j)+hape(i-1,j)+hape(i,j+1)+hape(i-1,j+1))
     end do
    end do
      
      hapCent    = 0.0
      hapScaleX  = 0.0
      hapScaleE  = 0.0
      hapCent = sqrt(hapcentx * hapcentx + hapcente * hapcente)    
      hapScaleX = sad * hapcentx/(1e-20 + hapCent) * rmask
      hapScaleE = sad * hapcente/(1e-20 + hapCent) * rmask

    ! Moved swim-all down below, Units here are fractions of 1.0
      GF%u_swim(spn,:,:) = hapScaleX ! hapcentx * swim_all  ! STILL in m/sec
      GF%v_swim(spn,:,:) = hapScaleE ! hapcente * swim_all ! STILL in m/sec

end do  ! spn loop

END SUBROUTINE feast_biol

!--------------------------------------------------------------------------------- 

SUBROUTINE fbiol_ind(NUM_ROMS_DEPTHS,rmask,ROMS_depth,ROMS_edges,&
        &ROMS_temp,ROMS_zoop,fal,fl,fsp,yday,LBi,UBi,LBj,UBj,in_spn,PPYY_out)

        USE mod_scalars
        USE mod_param
        
        IMPLICIT NONE 
        !TYPE (T_FEAST) :: GF
        integer, intent(in) :: LBi,UBi,LBj,UBj,NUM_ROMS_DEPTHS,in_spn
        real(r8), intent(in) :: ROMS_depth(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS)
        real(r8), intent(in) :: ROMS_edges(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS+1)
	      real(r8), intent(in) :: ROMS_temp(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS)
	      real(r8), intent(in) :: ROMS_zoop(LBi:UBi,LBj:UBj,NUM_ROMS_DEPTHS,NUM_PLANKTON)
	      real(r8), intent(in) :: rmask(LBi:,LBj:)

        real(r8), intent(inout) :: fal(LBi:UBi,LBj:UBj,nfvaral, &
                               &NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
        real(r8), intent(inout) :: fl(LBi:UBi,LBj:UBj,nfvarl, &
                               &NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS)
        real(r8), intent(inout) :: fsp(LBi:UBi,LBj:UBj,nfvar, &
                               &NUM_SIMPLE_SPECIES)
        real(r8),intent(in) :: yday
        real(r8), intent(out) :: PPYY_out(LBi:UBi,LBj:UBj)
        
      ! Internal indices
        integer  :: i,j,k,z,koff,fvaral,spal,lc,ac,fvarl,spl,fvar,sp
				integer  :: wrap_outflag,ictr,jctr
				
  ! by-layer cell properties
    real(r8) :: mixed
    real(r8), dimension(LAYERS)              :: in_T, in_D      
	  real(r8), dimension(LAYERS,NUM_PLANKTON) :: zoop   
	  real(r8), dimension(LAYERS,NUM_PLANKTON) :: zooptot   	  
    integer, dimension(NUM_ROMS_DEPTHS)      :: ROMS_layer  
    real(r8), dimension(NUM_ROMS_DEPTHS)     :: ROMS_layer_prop 
    
! Run_Step Declarations
  integer :: a,py,lpy,d,ip,pind,age
  real(r8), dimension(LAYERS) :: depth_pref 

              ! Equation 1
                real(r8) TT
				      ! Equation 2 and 2a - weight and biomass
                real(r8) B !,W
              ! Equation 6 - available prey and preference
              ! Preference from equation 7 computed elsewhere and saved
								real(r8) PPYY, zoop_pref !, Smax, Ratio
							! Equations 8 and 10 - diet composition and prey energy density	 
                real(r8) EDpy
                 ! Equation 9 - max. consumption as function of l, temperature
                real(r8) :: Cmax_stom, Cmax_bioe    !double Cmax;
							! Equation 5 // simple for now!
                !real(r8) ACT, Eact 
								!real(r8) :: Ar_W_Br, CC_search, CC_resp, Vel_ideal, Vel_actual								
          			real(r8) ::	pj, ASSIM, CperV, Vc, Fc_T, Cmax, Avel, Jmax
								real(r8) ::	Ar_W_Br, handle, Vr, Fr_T, Rfix, Rvel
								real(r16) :: cubeterm, cA2, cA1, cA0, cQQ, cRR, cDD, cSS, cTT, Vel_ideal
				        real(r8) ::  Cj_assim, Cpy, Cj, RvelV2, Rj, ACT_plus, SFU
							
							! Equations 3, 4, 10	  
                !real(r8) :: TotCons(TOT_FEAST+NUM_PLANKTON)
								real(r8) :: bioen_P, handle_P
                real(r8), dimension(TOT_LINKS) ::  tCpy, tED, tCj, tPP, tHP
                real(r8), dimension(TOT_LINKS,TOT_FEAST+NUM_PLANKTON) :: tTotCons
                ! For forcing vonB growth
                  real(r8) vonB_H, vonB_D, vonB_K, vonB_delW
                ! Equations 11, 12, 13
                  !real(r8) Cj, RFU
                  !real(r8) V, f_T
					      ! SPECIES AND LENGTH Loop II
                  real(r8) recCF												 								               
                  real(r8) Lzero, Lpow, dLen, Ltwo
                  real(r8), dimension(TOT_LINKS) :: delG, delZ, delCAL, delCF, delM, delRec, delSpawn
                  real(r8) propChange, propCAL, propCF, propL
									!real(r8) :: G_len(NUM_AGES), G_lenMinus1(NUM_AGES)
                ! recruitment
                  real(r8) divr, zers
                                              
!  These are values that are needed to pass between loops below.  
! Values needed only within loops are declared within those loops.
   real(r8), dimension(TOT_LINKS) :: mort 
   real(r8), dimension(TOT_LINKS) :: delW, delJ, growJ, spawnJ, totJ 
	 real(r8), dimension(TOT_LINKS) :: fishDC     
   real(r8), dimension(LAYERS,NUM_PLANKTON)    :: ZoopDC
   real(r8), dimension(LAYERS,NUM_PLANKTON)    :: ZoopMort
   real(r8), dimension(LAYERS,NUM_PLANKTON)    :: zoop_WW_eaten
   !real(r8), dimension(LAYERS,NUM_PLANKTON)    :: ZoopMortRate
   real(r8) :: ZoopSurviveRate
   real(r8), dimension(NUM_SIMPLE_SPECIES)     :: SimpDC,SimpleMort,simple_WW_eaten,SimpleMortRate
   real(r8), dimension(NUM_SIMPLE_SPECIES)     :: simple_Bio
   real(r8) :: Ntmp(1:MAX_LCLASS,1:NUM_AGES)
   real(r8) :: N_at_len(MAX_LCLASS)
! Happiness functions
    real(r8) :: happy(LBi:Ubi,LBj:UBj)
    real(r8) :: sad(LBi:Ubi,LBj:UBj)
    real(r8) :: hapXdn(LBi:Ubi,LBj:UBj)
    real(r8) :: hapEdn(LBi:Ubi,LBj:UBj)
    real(r8) :: landXdn(LBi:Ubi,LBj:UBj)
    real(r8) :: landEdn(LBi:Ubi,LBj:UBj)   
    real(r8) :: hapcentx(LBi:Ubi,LBj:UBj)
    real(r8) :: hapcente(LBi:Ubi,LBj:UBj)
    real(r8) :: hapCent(LBi:Ubi,LBj:UBj)
    real(r8) :: hapScaleX(LBi:Ubi,LBj:UBj)
    real(r8) :: hapScaleE(LBi:Ubi,LBj:UBj) 
    real(r8) :: hapx(LBi:Ubi,LBj:UBj) 
    real(r8) :: hape(LBi:Ubi,LBj:UBj) 
   integer NN, NS, gr, spn, spy, typeflag, MM,SS,YY,cc, ibm
   integer, dimension(NUM_GEARS*TOT_LINKS) :: gr_sp, gr_gr, gr_link
   integer, dimension(TOT_LINKS) :: murder
   real(r8), dimension(TOT_LINKS) :: N_all,CF_all,W_all,CAL_all  
   real(r8), dimension(NUM_GEARS,TOT_FEAST) :: CPUF, F_required, F_actual
   real real_day
					  
   real_day = yday

wrap_outflag = 0

NN   = ALL_LINKS
!GF%ghfal = 0.0
!GF%ghfl  = 0.0
!GF%ghfsp = 0.0
!GF%ozm   = 1.0

!GF%eat_all   = 0.0
!GF%eat_count = 0.0
!GF%mort_all  = 0.0
!GF%eat_plank = 0.0
!GF%eat_simp  = 0.0
!GF%gCmax=0.0; GF%gRfix=0.0; GF%gVel_ideal=0.0; GF%gRvelV2=0.0; GF%gACT=0.0; GF%gMzero=0.0
!GF%gPPYY=0.0; GF%gCperV=0.0; GF%gCpy=0.0; GF%gJmax=0.0; GF%gCj=0.0; GF%gRj=0.0; GF%gSFU=0.0
!GF%gEDpy=0.0; GF%gMtwo=0.0
  ictr=LBi+(UBi-LBi)/2
  jctr=LBj+(UBj-LBj)/2

 cc=0
! --------------------BEGIN MAIN I,J loop (lasts until end of subroutine)
PPYY_out = 0.0
 
DO i=LBi,UBi
 DO j=LBj,UBj    
  if (rmask(i,j)>0.1) then

  cc=cc+1
! Subroutine to set layers and set temperature and zooplankton concentrations!
! by layers
  call set_mixed(LBi,UBi,LBj,UBj,i,j, &
	     & NUM_ROMS_DEPTHS,ROMS_depth,ROMS_edges,ROMS_temp,ROMS_zoop, &
	     & mixed,in_T,in_D,zoop,zooptot,ROMS_layer,ROMS_layer_prop)
     
	! ZERO THE LOCAL ARRAYS
	     murder = 0				
			 mort = 0.0
			 zoopMort = 0.0
			 zoop_WW_eaten = 0.0
			 simpleMort = 0.0
			 simple_WW_eaten = 0.0
			 N_at_len = 0.0
       !out_catch = 0.0 
       delW      = 0.0
       CPUF = 0.0
       MM = 0       
       simple_bio = 0.0
       
  ! LOAD LOCAL FISH CONDITIONS INTO NN ARRAY
  do sp=1,NUM_SIMPLE_SPECIES
     simple_bio(sp) = fsp(i,j,1,sp)
  end do
  
  do spn=1,NN
     sp = sp_all(spn)
     lc = lc_all(spn)
     ac = age_all(spn)
          
     if (sp>NUM_AGED_SPECIES) then
        N_all(spn)   = fl(i,j,1,sp - NUM_AGED_SPECIES,lc)
        CF_all(spn)  = fl(i,j,2,sp - NUM_AGED_SPECIES,lc)
        CAL_all(spn) = fl(i,j,3,sp - NUM_AGED_SPECIES,lc)
		 else
        N_all(spn)   = fal(i,j,1,sp,lc,ac)
			  CF_all(spn)  = fal(i,j,2,sp,lc,ac)
        CAL_all(spn) = fal(i,j,3,sp,lc,ac)
		 end if
		 if (CF_all(spn)<F_EPSILON) then
		    CF_all(spn) = 1.0
		 end if
     W_all(spn)  = all_AL_LL_BL(spn) * CF_all(spn)
		 totJ(spn)   = W_all(spn) * CAL_all(spn)		 
     
  !  2015 INDIVIDUAL FUNCTION GEARS DELETED
	end do
 	
  !  2015 INDIVIDUAL FUNCTION GEARS DELETED
  
  !do SS = 1,NN    ! RUN THROUGH IN REVERSE LENGTH ORDER!
     spn = in_spn
     SS = rev_order(SPN)
         
     if (N_all(spn)> -1.) then
	
		   sp = sp_all(spn)   
     ! Zero link-specific quantities
       !PredDC  = 0.0
       fishDC  = 0.0
			 zoopDC  = 0.0
			 simpDC  = 0.0 						        
			 !TotCons = 0.0
     ! Set depth zone to middle if water column structure is completely mixed.
! NO LAYERS				 if (mixed == 1) THEN
					depth_pref(1)=0.0
					depth_pref(2)=1.0
					depth_pref(3)=0.0

     ! Equation 1
       TT = 0.0
       DO d=1,LAYERS
          TT = TT + depth_pref(d) * in_T(d)
       END DO
              
		! Equation 2 and 2a - weight and biomass
    	B = W_all(spn) * N_all(spn) 
                 
    ! Equation 6 - available prey and preference
    ! Preference from equation 7 computed elsewhere and saved
		! 10/28/11 Change to relative preference
			PPYY = 0.0
		! Non-plankton species
			do YY = SS,NN ! don't have to do this for species longer than this one
				 spy = all_order(YY) !do spy = 1,NN
				 py  = sp_all(spy)
			 ! pre-calculated and stored gamma functions for preferenc
			 !
			 	if (N_all(spy).gt.F_EPSILON) then	
				   fishDC(spy) = N_all(spy) * W_all(spy) * all_preference(spn,spy)
				   PPYY = PPYY + fishDC(spy)
				else
				   fishDC(spy) = 0.0
				end if 
			end do		 					    
	 ! Plankton species
		 DO py = 1,NUM_PLANKTON 
			DO d = 1,LAYERS 
			  if (zooptot(d,py).gt.F_EPSILON) then
				   ZoopDC(d,py) = zooptot(d,py) * depth_pref(d) * zoop_preference(spn,py) 
				   PPYY = PPYY + ZoopDC(d,py)          
				 else
				   ZoopDC(d,py) = 0.0
				 end if                
			END DO
		 END DO								    
	 ! Simple species
		 DO py = 1,NUM_SIMPLE_SPECIES
		   if (simple_Bio(py).gt.F_EPSILON) then
				SimpDC(py) = simple_Bio(py) * simple_preference(spn,py)
				PPYY = PPYY + SimpDC(py)
			 else
			  SimpDC(py) = 0.0
			 end if                         
		 END DO
  
	 ! Equations 8 and 10 - diet composition and prey energy density	 
     !if (ibm>0) GF%PPYY_ibm(i,j,ibm) = PPYY
		 EDpy=0.0
     if (PPYY<F_EPSILON) THEN
				EDpy=1.0    
     else            
			  do YY = SS,NN ! don't have to do this for species longer than this one
					 spy         = all_order(YY) !do spy = 1,NN
					 py          = sp_all(spy)
					 fishDC(spy) = fishDC(spy) / PPYY   
           EDpy        = EDpy + CAL_all(spy) * fishDC(spy)   
				end do   
				DO py = 1,NUM_PLANKTON    
				 DO d = 1,LAYERS   
					 ZoopDC(d,py) = ZoopDC(d,py) / PPYY  
					 EDpy         = EDpy + ZoopDC(d,py) * fsh_Zoop_JperG(py)  
				 END DO  
				END DO		
				DO py = 1,NUM_SIMPLE_SPECIES    
					 SimpDC(py) = SimpDC(py) / PPYY
					 EDpy       = EDpy + SimpDC(py) * fsh_Simple_JperG(py)  
				END DO		
		 end if

	! CUBIC variable activity model (Feb 1 2012)	

				pj       = PPYY * EDpy
				ASSIM    = 1.0 - fsh_P_j(sp)
        !A_W_B = A_L * (L_all(spn) ^ B_L)
        !Wt =A_W_B*CF_all(spn)
        !totJ     = W_all(spn) * CAL_all(spn)
        if (TT.ge.fsh_R_TM(sp)) then
           murder(spn) = 1
           Cpy  = 0.0
			     Cj   = 0.0
			     SFU  = 0.0
           Rj   = 0.0
           Cmax = -2.0
        else if (TT.ge.fsh_C_TM(sp)) then   
           Cpy  = 0.0
			     Cj   = 0.0
			     SFU  = 0.0
			     Cmax = -1.0
           Vr       = fsh_RTMoverRTM_RT0(sp) - TT / fsh_RTM_RT0(sp)
           Fr_T     = (Vr ** fsh_X_R(sp) ) * exp(fsh_X_R(sp) * (1.0 - Vr))
           Ar_W_Br  = fsh_A_R(sp) * (W_all(spn) ** fsh_B_R(sp) ) * 13560.0
					 Rj       = Ar_W_Br * Fr_T
			 else 
        CperV    = PPYY * all_Aenc_LL_Benc(spn) !*EncRate*(L_all(spn)^Enc_Lpow)        
        !Vc       = (CTMj-TT)/(CTMj-CTOj)
         Vc       = fsh_CTMoverCTM_CT0(sp) - TT / fsh_CTM_CT0(sp)
				 Fc_T     = (Vc ** fsh_X_C(sp) ) * exp(fsh_X_C(sp) * (1.0 - Vc))
         Vr       = fsh_RTMoverRTM_RT0(sp) - TT / fsh_RTM_RT0(sp)
         Fr_T     = (Vr ** fsh_X_R(sp) ) * exp(fsh_X_R(sp) * (1.0 - Vr))
         Cmax     = Fc_T * all_AS_LL_BS(spn) !cmax_a * L_all(spn)^cmax_b
         Avel     = EDpy * ASSIM * CperV
         Ar_W_Br  = fsh_A_R(sp) * (W_all(spn) ** fsh_B_R(sp) ) * 13560.0				
				 Rfix     = Ar_W_Br * Fr_T      
         Rvel     = all_Ar_v_LL(spn) * Rfix  !A_v =Av_min+B_Lzero/(1+(B_Lzero/B_Lone-1)*L_all(spn)^B_Lpow)        
				 if (Avel < 1.0) then
				   Vel_ideal = 0.0
				 else
				  Jmax     = EDpy * ASSIM * Cmax
          handle   = 1.0/Jmax
          cubeterm = 2.0*Rvel*handle*handle*Avel*Avel
          cA2      = 4.0*Rvel*handle*Avel/cubeterm
          cA1      = 2.0*Rvel/cubeterm
          cA0      = -Avel/cubeterm
          cQQ       = (3.0*cA1 - cA2*cA2)  / 9.0
          cRR       = (9.0*cA2*cA1 - 27.0*CA0- 2.0*cA2*cA2*cA2) / 54.0
          cDD       = cQQ*cQQ*cQQ + cRR*cRR
          cSS       = (cRR+(SQRT(cDD))) ** (0.3333333333333333333333)!^(1/3)
          cTT       = (cRR-(SQRT(cDD))) ** (0.3333333333333333333333)!^(1/3)
          Vel_ideal = -(cA2/3.0) + (cSS+cTT)
         end if				  
				Cj_assim  = EDpy * ASSIM * CperV * Vel_ideal/ &
				          & (1.0 + handle * EDpy * ASSIM * CperV * Vel_ideal)
        Cpy       = Cj_assim / (EDpy*ASSIM)
        Cj        = Cpy * EDpy
        RvelV2    = Rvel*Vel_ideal*Vel_ideal
        Rj        = Rfix + RvelV2
        ACT_plus  = RvelV2/Rfix
        SFU       = fsh_P_j(sp) * Cj
			end if
			 
      ! recruitment
        delSpawn(spn) = 0.0
        if (yday.ge.fsh_sp_sday(sp) .and. yday.le.fsh_sp_eday(sp) ) then
             divr = (fsh_sp_eday(sp) - fsh_sp_sday(sp) + 1)
             delSpawn(spn)   = W_all(spn) * all_mat_rec(spn) / divr
        end if
        
      ! Summing the joules  
        growJ(spn)  = Cj - Rj - SFU
        spawnJ(spn) = delSpawn(spn) * CAL_all(spn)
      							
       ! Equation 14 predation mortality
				 do YY = SS,NN ! don't have to do this for species longer than this one
						spy=all_order(YY) !do spy = 1,NN                             
            mort(spy) = mort(spy) + N_all(spn) * Cpy * fishDC(spy)/W_all(spy)
         END DO    
								             
        ! This conversion is for trawl survey data per hectare
          ! Units:  DC * g/pred * pred/10,000m^2(1 hectare) * 1/depth(m) * 1000mg/g  = mg/m^3
        ! This one is for pred/m^2
          ! Units:  DC * g/pred * pred/m^2(1 hectare) * 1/depth(m) * 1000mg/g  = mg/m^3
          DO py = 1,NUM_PLANKTON    !for (py=0; py<NUM_PLANKTON; py++){
             DO d = 1, LAYERS
                 if (in_D(d)>0) THEN
                     zoop_WW_eaten(d,py) = zoop_WW_eaten(d,py) + ZoopDC(d,py) * Cpy * N_all(spn) 									 
								 else
               	     zoop_WW_eaten(d,py) = 0.0
								 end if												
             END DO								
         END DO				

        ! This loop is for the simple (not-lenghted) FEAST species
          DO py = 1,NUM_SIMPLE_SPECIES    !for (py=0; py<NUM_PLANKTON; py++){
             simple_WW_eaten(py) = simple_WW_eaten(py) + SimpDC(py) * Cpy * N_all(spn)				
         END DO	
     	    
							 tCpy(spn) = Cpy
							 tED(spn)  = EDpy
							 tCj(spn)  = Cj
							 tPP(spn)  = bioen_P
							 tHP(spn)  = handle_P
					 		 
        end if  !end of if N>0 loop
      
!	END DO!} //end of first species spn loop -------------------------------------

    end if ! END OF RMASK IF
    
    PPYY_out(i,j) = PPYY
    
  END DO ! main j loop
END DO  ! main i loop

END SUBROUTINE fbiol_ind


#endif 

END MODULE mod_feast
