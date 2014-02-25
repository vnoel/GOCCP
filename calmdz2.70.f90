!****************************************************************************!
!                                                                            !
!                            GOCCP2.55 same as 2.50 with false liq and ice   !
!                            files with occ with good discrimination threshold
!
!                    Gregory CESANA / LMD / IPSL                             !
!                      last update 27/05/2012                                !
!                                                                            !
!                                                                            !
!  1) Purpose : compute instantaneous and daily mean profiles of Scattering  !
!               Ration (SR), Color Ratio (CR) & Depolarization Ratio (DR)    !
!               over a model grid.                                           !
!                                                                            !
!               Compute various daily/monthly cloudiness over a model grid : !
!                  - Map of Low Mid High cloud cover                         !
!                  - 3D Cloud Fraction                                       !
!                  - 3D Cloud Phase                                          !
!                  - SR Histograms                                           !
!                                                                            !
!  2) Input   : SDS & META variables from CALIPSO hdf level1 Data files.     !
!                                                                            !
!  3) Output  : - 3D_CloudFraction :                                         !
!                    clcalipso(lon,lat,alt,time)                             !
!                    clrcalipso(lon,lat,alt,time)                            !
!                    uncalipso(lon,lat,alt,time)                             !
!                                                                            !
!               - Map of Low Mid High :                                      !
!                    cllcalipso(lon,lat,time)                                !
!                    clmcalipso(lon,lat,time)                                !
!                    clhcalipso(lon,lat,time)                                !
!                    cltcalipso(lon,lat,time)                                !
!                    clccalipso(lon,lat,time)                                !
!                                                                            !
!               - SR_Histo :                                                 !
!                    cfad_lidarsr532_Occ(lon,lat,alt,box,time)               !
!                    cfad_lidarsr532_Occ2(lon,lat,alt,box2,time)             !
!                                                                            !
!               - 3D_CloudPhase :                                            !
!                    ice_cloud(lon,lat,alt,time)                             !
!                    water_cloud(lon,lat,alt,time)                           !
!                                                                            !
!               - instant_SR_CR_DR                                           !
!                    longitude(it)                                           !
!                    latitude(it)                                            !
!                    altitude(it)                                            !
!                    time(it)                                                !
!                    SE(it)                                                  !
!                    instant_SR(it,alt)                                      !
!                    instant_CR(it,alt)                                      !
!                    instant_DR(it,alt)                                      !
!                                                                            !
!  4) Grid :                                                                 !
!               -CFMIP2   : 2° x 2° x 40levels from 0 to 19.2km              !
!                           (180,90,40)                                      !
!               -CFMIP1   : 1° x 1° x 40levels from 0 to 19.2km              !
!                           (360,180,40)                                     !
!               -CFMIP2.5 : 2.5° x 2.5° x 40levels from 0 to 19.2km          !
!                           (144,72,40)                                      !
!               -CFMIP    : 3.75° x 2.5° x 40levels from 0 to 19.2km         !
!                           (96,72,40)                                       !
!               -LMDZ     : 3.75° x 2.5° x 19levels from 0 to 40.4823km      !
!                           (96,72,19)                                       !
!               -LMDZ40   : 3.75° x 2.5° x 40levels from 0 to 40.4823km      !
!                           (96,72,40)                                       !
!               -NASA     : 5° x 5° x 40levels each 480m from 0 to 19.2km    !
!                           (73,37,41)                                       !
!                                                                            !
!  5) Compilation : use the makefile "makefiles.sh"                          !
!      makefiles.sh $1                                                       !
!      ifort $1.f90 -I/usr/include/hdf -L/usr/lib64/hdf -lmfhdf -ldf -ljpeg  !
!                   -lz -I/opt/netcdf-3.6.0-p1-ifort-64/include              !
!                   -L/opt/netcdf-3.6.0-p1-ifort-64/lib  -lnetcdf -o $1.e    !
!                                                                            !
!  6) Last updates & bug fix :                                               !
!   - 02/02/08 altitude files changed                                        !
!   - 05/02/08 add CFMIP horinzontal grid                                    !
!   - 06/02/08 add NASA horizontal grid                                      !
!   - 21/03/08 add PSC filter                                                !
!   - 29/03/08 memory optimization                                           !
!   - 12/04/08 bug in atb_mol fixed                                          !
!   - 15/04/08 bug in filtre_2lvl fixed                                      !
!   - 21/05/08 in progress : pressure mode                                   !
!   - 21/06/08 documentation                                                 !
!   - 21/06/08 maj subgrid3                                                  !
!   - 14/09/08 change in diagSR calculation                                  !
!   - 15/09/08 change in diagSR boxes from -887 to 100000                    !
!   - 16/09/08 add monthdiagSR_Occ & monthdiagSR_Frac in the diagSR files    !
!   - 15/10/08 change SrSeuilCloud from 3 to 5                               !
!   - 17/10/08 add negdiagSR_Occ + new diag boxes :                          !
!             -887,0,0.01,1.2,3,5,7,10,15,20,25,30,40,50,60,80,10e05         !
!   - 29/01/09 add Seuil_rap 0.25e-28 0.97e-28                               !
!   - 06/03/09 longitude latitude changées                                   !
!   - 23/03/09 longitude latitude altitude fixed                             !
!   - 23/03/09 add LMDZ40 lvl                                                !
!   - 23/03/09 add lon,lat,alt middle in the outputs                         !
!   - 23/03/09 attribut variable s name changed in output files              !
!   - 23/03/09 log version of the program added                              !
!   - 24/03/09 add cloud filter over 21km to LMDZ40                          !
!   - 01/04/09 bugfix topmidl = 15 for CFMIP GRID                            !
!   - 01/04/09 exclude sat fraction of the cloud fraction calcul             !
!   - 01/04/09 exclude fully attenuated layer in the calcul of low mid high  !
!              cloud type.                                                   !
!   - 19/04/09 change file's name output                                     !
!   - 22/04/09 add indice to calculate nb of retrievals and number of        !
!              observations in maplowmidhigh files.                          !
!   - 28/05/09 add CFMIP80 grid                                              !
!   - 01/06/09 bugfix on SR_hist calculation (initialization monthdiagSR1)   !
!   - 09/06/09 add calculation of Color Ratio                                !
!   - july 09  change in the lat/lon files (from negative value to positive) !
!   - july 09  add ntot (total number of point used), and nret, nretl, nretm,!
!              nreth (number of retrieval tot low mid high)                  !
!   - 10/09/09 bug fix on low with the selection of profil non fully         !
!              attenuated, (treshold fixed)                                  !
!   - 10/09/09 add CFMIP 2deg grid                                           !
!   - 29/10/09 add WRF grid                                                  !
!   - 29/10/09 add DEPOL CR                                                  !
!   - 29/10/09 add new srbox in depol                                        !
!   - 29/10/09 add parameters file to lighten the code                       !
!   - 29/10/09 add time controle                                             !
!   - change in "select type of hdf file" to match with version 3 of hdf file!
!   - 13/01/10 add empty test on input file list and hdf input files         !
!   - 25/01/10 add bugfix in colclear calcul                                 !
!   - 25/01/10 add srbox_bound                                               !
!   - 26/01/10 change in topmid lvl                                          !
!   - 26/01/10 bugfix in sat layer low mid high lvl(8-10% less middle cloud) !
!   - 23/02/10 depol saved in .ascii file                                    !
!   - 23/02/10 optimization of several loop and parts of the code such as    !
!              use of exit instead of goto...                                !
!   - 23/02/10 change real*8 and real*4 var to integer or real               !
!   - 23/02/10 change in surf_detect2 -888 instead of -9999                  !
!   - 23/02/10 change in zero_detect -777 instead of -9999                   !
!   - 23/02/10 change in SR_CR_DEPOL_mean add -777 and -9999                 !
!   - 23/02/10 the same for filtre_2lvl                                      !
!   - 23/02/10 SE_detect 1 2 change in the selection of first level + exit   !
!              instead of goto                                               !
!   - 23/02/10 add rej_fraction in fraction_subgrid2                         !
!   - 23/02/10 add indice for nan value in the SR_histo                      !
!   - 23/02/10 change length of srmod by adding 2 boxes -888 and -777        !
!   - 23/02/10 add -777 -888 in instant_SR routine                           !
!   - 02/03/10 modif in the SR_histo recording (suppression of monthdiagSR1  !
!              and monthdiagSR.                                              !
!   - 03/03/10 bugfix in the calculation of the 3DCF indice (add rejfraction)!
!   - 19/03/10 change in the SR_histo files, add of cfad2 with the three     !
!              firts boxes of old cfad. Modif of old cfad from 18 box to 15  !
!   - 20/04/10 add switch in fraction_subgrid3                               !
!   - 21/04/10 bugfix in the routine surf_detect2 if alt=SE then -888        !
!   - 28/06/10 bugfix srbox_bound                                            !
!   - 12/07/10 add instant_SR2 routine = instant_SR corrected by delta ATB   !
!   - 12/07/10 add instant_SR_DEPOL routine                                  !
!   - 26/07/10 add CFMIP16 and CFMIP32 vertical grid 160 & 320 levels        !
!   - 02/09/10 Change in mid top level from 7.2 to 6.5km                     !
!   - 15/10/10 Change SR_DEPOL routine into SR_CR_DEPOL : Calculation of SR, !
!              CR & DEPOL instantaneous variables                            !
!   - 15/12/10 Add alt_mid alt_bound in instant_SR files                     !
!              Bugfix in instant_SR file file date                           !
!   - 16/12/10 Add fraction_subgrid_8km = improve the SNR for day time data  !
!   - 17/12/10 Add new variables for the GEWEX cloud assessment              !
!                 - Temperature of cloud (tot low mid high)                  !
!                 - Averaged altitude of cloud                               !
!                 - Histograms of 2D var (temp, alt, cl tot/low/mid/high)    !
!              ==> add create_mapnc2 map_recvar2nc3 map_recvar2nc4 routines  !
!   - 20/12/10 Add 3D Cloud Phase file with create_depolnc3D &               !
!              depol_recvar2nc routines                                      !
!   - 20/12/10 Add vertical_mean_hori routine to horizontally average the    !
!              variables                                                     !
!   - 05/01/11 Change in references of all files                             !
!   - 05/01/11 Add alt_mid alt_bound in all files                            !
!   - 05/01/11 Add Color Ratio variables in instant files                    !
!   - 05/01/11 Bugfix in fraction_subgrid_8km in the upper/lower bound of the!
!              flag calculation                                              !
!   - 06/01/11 Unification of Gewex files in a single one                    !
!   - 06/01/11 Substitution of instant_SR files to instant_SR_CR_DR files    !
!   - 25/02/11 Bug in indphasemonth
!   - 02/03/11 add -777 value depending on delta to srmoy before recording
!              instantaneous SR files
!              ==> Major change in SR_histo files during DAY time
!   - 02/03/11 add -777 value to DEPOL and CR where srmoy=-777 (except for 
!              the delta case)
!   - 02/03/11 add 2dim srbox in SR_histo files
!   - 16/03/11 bugfix on the instant_SR_CR_DR filename
!   - 16/08/11 change in SR seuil cld in fraction subgrid 8km routine
!              from 13.3 to 15
!   - 17/08/11 add commande system to change instantSR name
!              add routine fraction_subgrid3_8km for 1-7 flag figs
!              add fraction_subgrid2_8km_delta
!   - 25/08/11 add computation of ATB ATBr ATBl & ATBmol after 
!              fraction_subgrid 
!   - 25/08/11 add cloud phase diag based on ATBr=f(ATB) relation
!   - 27/08/11 add routine SR_CR_DR_ATB_nc  Same routine as SR_CR_DR_2nc 
!              including the record of ATB ATBper ATBpar & ATBmol variables
!   - 31/08/11 calculation of CP3D in the same way as CF3D
!              add indphaseday indphasepermonth and filter under 8.16km of
!              cloud below cloud with SR>30
!   - 26/10/11 add new var: isccp low mid high for UNDEF/HO and dust cloud
!              add var 3D UNDEF/HO/DUST phase cloud
!   - 30/11/11 add non over lap mode if nol = 1 nol activated
!   - 07/12/11 add switch to select type of instant sr file
!              add calcul of seuilsnrlow and high for fraction_subgrid_8km
!              (day time) to be adapted to any kind of grid
!   - 11/12/11 new tempmod grid from -80 to 30 degree every 2degrees
!   - 23/01/12 correction ind routine interp: boundaries of the domain were
!              excluded with -gt instead of -ge
!              ===> change in molecular and temperature profiles
!              ===> statistics change over land
!   - 23/01/12 add nanlow when no low profiles go through the low layer
!              ===> nan value
!   - 24/01/12 correction of atb_mol_interp routine in boundaries of the
!              domain
!              add switch option in atb_mol_interp routine to interp 
!              temperature close to the ground
!              add temperature criteria for dust case between -45 & 45 lat
!              for cloud above the discrimination threshold and below 3km
!              add cloud with depol > 1 in undefined phase type
!   - 06/02/12 correction isccpdustday out of nphase loop
!              correction in dimension of isccpunday in nphase loop
!   - 06/02/12 change in temperature criteria for dust case between -45 & 45 
!              of latitude to all latitude cases -90 to 90 and for cloud above 
!              the discrimination threshold and below 5km instead of 3km
!   - 06/02/12 hocloud becomes cloud flagged as liquid with temperature below
!              -41°C
!   - 06/02/12 add instant_switch variable to active instant_SR file or not
!   - 06/02/12 add routine instant_phase() to record instantaneous phase in 
!              instant_SR_CR_DR files.
!              instant_phase is ranged from 1 to 5:
!               1=LIQ / 2=ICE / 3=UNDEF / 4=FakeICE / 5=FakeLIQ
!   - 08/02/12 add fraction phase (percentage of water/ice content) in 3D & 
!              Map files
!   - 14/02/12 correction of ground lidar echo until the thirth lvl above the 
!              ground in SE_alt_atb routine
!   - 21/03/12 add version in name files automatically
!              change in variable names and files name for output phase 
!              files
!              all unclassified phase variables (3D and MapLMH) are unified 
!              in one variable with 5 categories for: 
!              1=undefined phase e.g cloud below cloud with SR>30
!              2=fake ice e.g cloud flaged as ice but with TEMP>0°C
!              3=fake liq e.g cloud flaged as liq but with TEMP<-42°C
!              4=HO e.g horizontally oriented ice particle cloud
!              5=Noise e.g cloud with DEPOL>1 (physically impossible)
!              Correction in URL of the website in output files
!              Add cloud mask and cloud phase in instantaneous files
!   - 28/03/12 Add PSCs filter during Arctic winter +10% profile when 
!              PSCs occurs 
!              old version: for month 201001 2.6% of nan profiles
!              new version: for month 201001 0.8% of nan profiles            
!   - 29/03/12 Correction in zero_detect routine, substitution of -9999 value
!              by -777
!    changement dans interp_mol un ge en gt
!    where srmoy = -777 tempmoy /= -777
!    surf_detect lt = le
!    zero detect -777 = -9999
!    vertical_mean ajout de -777
!    atb_mol -777
!    fraction_subgrid3_8km  fraction_subgrid2_8km  add rejfrac qd deltATB<Seuil
!    virer srmoy = -777 a cause du deltatb
!    ajouter where rej_fraction pour srmoy
!    ajouter filtre profil -777 avant vertical mean  
!    tempmod35 (35 boite de 3deg de chaque cote de 0 degC)
!    boite de 2deg trop petite = bug CFice et liq stratifié
!    vertical_mean special temperature nvar=4
!    2.43 change in SE_alt_atb, search ground lvl until 30lidar lvl above the
!    Surface Elevation (instead of 6 lvl in previous version) 
!    First lvl below cloud with SR>30 = phase of cloud above
!    "false liq"=ice & "false ice"=liq in statistics  
!    tempmax=39 add -90 ===> -78  
!   - 10/05/12 correction bug on temperature calculation. No -9999 and -777
!              values allocated only -888 
!   - 27/05/12 change in phase parameter description in instant_phase routine
!   - 27/05/12 Name of variable change in phase file, PIC becomes RPIC
!   - 03/06/12 Add control of nan value of SE in SE_alt_atb
!   - 06/07/12 Correction of atb_mol routine: when one -777 or -9999 value is
!              detected in the profile, all the profile is set to -777
!              old version had problems with nan molecular. 
!   - 02/10/12 add good phase discrimination threshold
!   - 02/10/12 add CF3D temp & phase files with occ (to estimate number of 
!              cases using T criterion 
!   - 02/10/12 retour au mode SR30 sans mettre la phase identique au nuage du
!              dessus
!   - 02/10/12 correction des undef dans les fichiers map et CF3D en ajoutant
!              les 5 cas (unclass, false liq, false ice, HO, unphysical) au 
!              lieux de 3 avant
!   - 12/10/12 Nouveau seuil de discrimination phase basé sur les stats JFM
!              overlapped
!   - 01/11/12 Add top/base cloud height of high cloud
!                                                                            !
!****************************************************************************!

!************************* SUBROUTINE SCHEME ********************************!
!                                                                            !
! sdsread8            : x2 for time                                          !
! sdsread             : x8 for atb,atb1064,mol,perp,temp,lat,lon,SE          !
! metaread            : x2 for altl & altm                                   !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! interp              : x4 for Pressure and molecular                        !
! atb_mol             : x4 for normalized ratio & molecular calculation      !
! atb_mol_interp      : x2 for molecular extrapolation                       !
! SE_alt_chim         : x4 for adding the SE to atb & mol                    !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! vertical_mean       : x6 for the average of atb, mol, perp, temp           !
! zero_detect         : x8 for the detection of empty boxes                  !
! Surf_detect2        : x8 for adding the Surface Elevation                  !
! SE_km_2_pres2       : x2 for adding the Surface Elevation                  !
! SR_CR_DEPOL_mean    : x16 for the SR CR & DR calculation                   !
! filtre_2lvl         : x4 for the SR calculation with a filter up to 21km   !
!                          (over 21kilometers grid such as LMDZ/LMDZ40 grid) !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! SR_CR_DR_2nc        : x1 for the recording of instant_SR instant_CR and    !
!                       instant_DR variables in netcdf files                 !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! fraction_subgrid2_8km   : x1 for cloud,clear,uncertain flags calculation   !
!                                                                            !
!----------------------------------------------------------------------------!
!                                                                            !
! create_mapnc        : x1 for the creation of the MapLowMidHigh file        !
! create_mapnc2       : x1 for the creation of the GEWEX MapLowMidHigh file  !
! map_recvar2nc2      : x1 for the recording of map variables                !
! map_recvar2nc7      : x1 for the recording of GEWEX map variables          !
! create_diagnc       : x1 for the creation of the SR_histograms file        !
! diag_recvar2nc3     : x1 for the recording of the SR_histograms variables  !
! create_profnc       : x1 for the creation of the 3D_CloudFraction files    !
! prof_recvar2nc      : x1 for the recording of the 3D_CloudFraction var     !
! create_depolnc3D    : x1 for the creation of the 3D_CloudPhase files       !
! depol_recvar2nc     : x1 for the recording of the 3D_CloudPhase var        !
!                                                                            !
!----------------------------------------------------------------------------!

program calmdz

  implicit none


!****************************************************************************!
!*!!!!!!!!!!!! DEFINITIONS & DECLARATIONS OF VARIABLES !!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


!----------------------------------------------------------------------------!
! Date and time variable                                                     !
!----------------------------------------------------------------------------!
!               variables used in subroutine DATE_AND_TIME()                 !
!****************************************************************************!
!   value     : 8 dimension value which contains the year,month,day,hour,    !
!               minute, seconds and milliseconds of the real-time.           !
!   datenc    : date from the real-system clock, format: yyyymmdd.           !
!   timenc    : time from the real-system clock, format: hhmmss.sss.         !
!   zone      : represente the difference with respect to Coordinated        !
!               Universal Time (UTC), format: (+-)hhmm.                      !
!                                                                            !
!----------------------------------------------------------------------------!
!               variables used all along the program                         !
!****************************************************************************!
!   date      : first value of the variable time in the CALIPSO hdf file.    !
!   year      : year of the CALIPSO hdf file (derived from date var).        !
!   month     : month of the CALIPSO hdf file (derived from date var).       !
!   day       : day of the CALIPSO hd file (derived from date & month var).  !
!   jour      : day of the CALIPSO hd file (derived from date & month var).  !
!   resd      : number of days since 2000/01/01 for the trimonth period.     !
!   resh      : number of hour since 2000/01/01 for the trimonth period.     !
!                                                                            !
!----------------------------------------------------------------------------!
!               variables used in the function interdat()                    !
!****************************************************************************!
!   interdat  : routine which allow to perform the calculation of the number !
!               of days since the 2000/01/01 to the middle of the period     !
!               looked by the routine.                                       !
!   date1     : input of the routine interdat. It is the date of the midle   !
!               on the period looked by the routine, format: yyyymmddhh      !
!   date2     : input of the routine interdat. It is the date of the         !
!               first day since the year 2000, format: yyyymmddhh.           !
!                                                                            !
!----------------------------------------------------------------------------!
      integer,dimension(8)  ::  value
      character  ::  datenc*8, timenc*10, zone*5 
      character  ::  blanks*15,forme*15
      character(len=15),parameter ::  version="2.70"
      integer  ::  date, year, month, day, jour          
      integer  ::  date1, date2, resh                  
      integer  ::  interdat, ret1
      real  ::  resd, test    

!----------------------------------------------------------------------------!
! Calculation of time elpased                                                !
!----------------------------------------------------------------------------!
      integer  ::  t1, t2, ir
      real     ::  t_cpu, t_cpu_0, t_cpu_1, tempstot
 
!----------------------------------------------------------------------------!
! Name of keyboarding data and some char variables                           !
!----------------------------------------------------------------------------!
!                            input files                                     !
!****************************************************************************!
!   file2     : path and name of hdf CALIPSO input file                      !
!   file3     : path and name of the file which lists all the hdf CALIPSO    !
!               file name and path (=file2)                                  !
!                                                                            !
!----------------------------------------------------------------------------!
!                            output files                                    !
!****************************************************************************!
!   file4     : instant_SR_CR_DR output file name                            !
!   file5     : 3D_CloudFraction output file name                            !
!   file6     : MapLowMidHigh output file name                               !
!   file7     : SR_histo output file name                                    !
!   file10    : 3D_CloudFraction_Phase output file name                               !
!   file11    : MapLowMidHigh_Phase output file name                               !
!                                                                            !
!----------------------------------------------------------------------------!
!                             model grid                                     !
!****************************************************************************!
!   model     : name of model selected (at that time only lmdz is available) !
!   gcm       : name of the grid selected description in section 4)          !
!                                                                            !
!----------------------------------------------------------------------------!
!                           switch variables                                 !
!****************************************************************************!
!   sauve     : = "chim" "wrf" or "off" to record output files in ASCII      !
!   switch    : = "night" or "day" to select day or night mode               !
!   switch2   : = "sat" or "cloudy" to select sat or cloudy mode             !
!   alt_pres  : = "altitude" or "pressure" to select altitude or pressure    !
!               mode.                                                        !
!                                                                            !
!----------------------------------------------------------------------------!
!                          other char variables                              !
!****************************************************************************!
!   metal     : name of lidar altitude variable metadata over 583lvl         !
!   metam     : name of meteo altitude variable metadata over 33lvl          !
!   sds_varname : name of SDS variables from hdf CALIPSO files               !
!   datec     : date in characters converted from date                       !
!   numfichc  : number of hdf files read in characters, converted from       !
!               numfich                                                      !
!                                                                            !
!----------------------------------------------------------------------------!
      character  ::  file2*1024,  file3*132, file4*132 , file5*132   
      character  ::  file6*132, file7*132, file8*132, file9*132
      character  ::  file10*132, file11*132, file12*132,file13*132 
      character  ::  file66*132   
      character  ::  model*30, filetmp*1024, filetmp2*1024,instantname*1024
      character  ::  metal*30,metam*30
      character  ::  command*1024, command2*1024, command3*1024,command4*1024
      character  ::  sauve*30,gcm*8, idiagc*2, idepc*2
      character  ::  sds_varname*100
      character  ::  switch*5, alt_pres*8, switch2*6,instant_switch*8
      character  ::  datec*6, numfichc*4, yearc*4
      character  ::  grid*8, altfile*10, lonfile*18, latfile*17
  100 format(A200)
                 
                                       
!----------------------------------------------------------------------------!
! id used in netcdf output files with subroutine create_*nc()                !
!----------------------------------------------------------------------------!
!   dimidsp   : dimension id of the MaP3D variables recorded in the ncdf     !
!               files                                                        !
!   dimidsm   : dimension id of the MaPLowMidHigh variables recorded in the  !
!               ncdf files                                                   !
!   dimidsd   : dimension id of diagSR variable recorded in the ncdf file    !
!   nan       : nan is the NaN value in the ncdf files                       !
!                                                                            !
!----------------------------------------------------------------------------!
      real,parameter  ::   nan=-9999.
      real*4,parameter  ::  SeuilMol1km = 0.00015 , SeuilTemp1km = 6.5
      real*4  ::  seuilatb, deltatb  
      integer  ::  dimidsp(4),dimidsp2(5), dimidsm(3), dimidsm2(4),dimidphase(4),dimidsd(5),dimidsdb(5)
      integer  ::  dimidhist3(4), dimidsd2(4),dimidsd3(5),dimidsd4(6),dimidpha(6)
      integer  ::  dimidhist(4),dimidhist2(4)

!----------------------------------------------------------------------------!
! Output variable dimensions                                                 !
!----------------------------------------------------------------------------!
!   latmax    : number of the latitude boxes, according to the model chosen  !
!               CFMIP=181, LMDZ=72, NASA=37                                  !
!   lonmax    : number of the longitude boxes, CFMIP=361, LMDZ=96, NASA=73   !
!   altmax    : number of the altitude boxes, CFMIP=41, LMDZ=19, NASA=41     !
!   diagmax   : number of the diagSR boxes, always equal to 15               !
!   daymax    : maximum number of day in a month (for the daily var)         !
!                                                                            !
!----------------------------------------------------------------------------!
      integer  ::  latmax, lonmax, altmax, altmax2
      integer, parameter  ::  diagmax = 19 , diagmax2 = 11, depolmax = 21
      integer, parameter  ::  daymax = 31, pr2max=301, permax=241, tempmax=39
      integer  ::  toplowl,topmidl,tophighl,toplvlsat1,toplvlsat2
      integer  ::  altstart,altend,nol = 0
                  
                     
!----------------------------------------------------------------------------!
! Loop variables                                                             !
!----------------------------------------------------------------------------!
!                        obs & model variables                               !
!****************************************************************************!
!   i         : loop index of the number of profil in the hdf file read      !
!   iz        : loop index of the altitude boxes                             !
!   altitude  : number of level of the lidar variables                       !
!                                                                            !
!----------------------------------------------------------------------------!
!                            obs variables                                   !
!****************************************************************************!
!   ilid      : loop index of the lidar altitude (583)                       !
!   altitude2 : number of level of the meteo variables                       !
!                                                                            !
!----------------------------------------------------------------------------!
!                        CHIMERE & WRF variables                             !
!****************************************************************************!
!   j         : loop index of the latitude boxes                             !
!   k         : loop index of the longitude boxes                            !
!                                                                            !
!----------------------------------------------------------------------------!
!                    model variables (daily/monthly)                         !
!****************************************************************************!
!   ilon      : loop index of the longitude boxes in the daily/monthly var   !
!   ilat      : loop index of the latitude boxes in the daily/monthly var    !
!   idiag     : loop index of the diagSR boxes in the daily/monthly var      !
!                                                                            !
!----------------------------------------------------------------------------!
!                         Unspecific variables                               !
!****************************************************************************!
!   numfich   : number of hdf file read, = number of line in file3           !
!   comptpf   : number of boxes crossed by the satellite                     !
!   box       : number total of boxes in the WRF or CHIMERE grid             !
!                                                                            !
!----------------------------------------------------------------------------!
      integer  ::  iz, i, j, k, ilat, ilon, ialt, ilid, idiag, idep, ipr2,iperp,itemp  
      integer(kind=2),parameter :: altitude2 = 33  
      integer(kind=2),parameter :: altitude = 583    
      integer  ::  numfich, comptpf,  box , nphase       

   
!----------------------------------------------------------------------------!
! Status variables                                                           !
!----------------------------------------------------------------------------!
!   err       : return 0 if the file could be opened                         !
!   OK_buffer : return 0 if the allocation is OK                             !
!                                                                            !
!----------------------------------------------------------------------------!
      integer  ::  err , OK_buffer                     
                     

!----------------------------------------------------------------------------!
! Input variable dimensions                                                  !
!----------------------------------------------------------------------------!
!   it        : number of profil (=nprof) in the hdf file read (about 60000) !
!                                                                            !
!----------------------------------------------------------------------------!
      integer :: it
      integer*4,dimension(:,:),allocatable :: indretmean, indtotmean, indtot
      integer*4,dimension(:,:),allocatable :: indretlowmean,indretmidmean
      integer*4,dimension(:,:),allocatable :: indrethighmean 
      integer*4,dimension(:),allocatable :: indret, indretlow, indretmid
      integer*4,dimension(:),allocatable :: indrethigh

                     
!----------------------------------------------------------------------------!
! SDS variables from CALIPSO hdf files, extracted by routine sdsread (1bit) &!
! sdsread8 (8bit)                                                            !
!----------------------------------------------------------------------------!
!   atb       : Total_Attenuated_Backscatter_532 per meter per steradian !
!               has dimension (altitude,nprof)                               !
!   mol       : Molecular_Number_Density in count numerical CN has dimension !
!               (altitude2,nprof)                                            !
!   pres      : Pressure in hPa has dimension (altitude2)                    !
!   temps     : Profil_UTC_time has format yymmdd.ffffffff, has dimension    !
!               (nprof)                                                      !
!   vartmp8   : Temporal var for the time reading, has dimension (nprof)     !
!   vartmp    : Temporal var for the lat,lon,SE reading, has dimension       !
!               (nprof)                                                      !
!   temps2    : Time converted in UTC fractionned hour of a day (ex : 1.5 =  !
!               1h30 am), has dimension (nprof)                              !
!   lat       : Latitude in degrees, has dimension (nprof)                   !
!   lon       : Longitude in degrees, has dimension (nprof)                  !
!   SE        : Surface_Elevation in kilometers, has dimension (nprof)       !
!                                                                            !
!----------------------------------------------------------------------------!
      real,dimension(:,:),allocatable  ::  atb, mol,  pres, atb2, perp, temp
      real*8,dimension(:),allocatable  ::  temps, temps2,temps2wrf
      real*8,dimension(:,:),allocatable  :: vartmp8
      real,dimension(:,:),allocatable  :: vartmp
      real,dimension(:),allocatable  ::  lat, lon
      real,dimension(:),allocatable  ::  SE,SEwrf
      real,dimension(:),allocatable  ::  latwrf,lonwrf
                    
!----------------------------------------------------------------------------!
! SDS variables interpolated verticaly from CALIPSO variables                !
!----------------------------------------------------------------------------!
!   mol2      : Molecular interpolated from 33lvl to 583lvl in CN            !
!               dim=(altitude,nprof)                                         !
!   pres2     : Pressure interpolated from 33lvl to 583lvl in hPa            !
!               dim=(altitude,nprof)                                         !
!   mol3      : Molecular converted from mol2(CN) in km-1 sr-1               !
!               dim=(altitude,nprof)                                         !
!   SEp       : Surface_Elevation in hPa derived from SE(km), dim=(nprof)    !
!                                                                            !
!----------------------------------------------------------------------------!
      real,dimension(:,:),allocatable  ::  mol2, pres2, mol3, mol4, temp2 
      real,dimension(:),allocatable  ::  SEp
                  
                   
!----------------------------------------------------------------------------!
! Output & Grid variables                                                    !
!----------------------------------------------------------------------------!
!   latmod    : values of model latitude from 90 to -90                      !
!               dimension : CFMIP=181, LMDZ=72, NASA=37                      !
!   lonmod    : values of model longitude from 180 to -180                   !
!               dimension : CFMIP=361, LMDZ=96, NASA=73                      !
!   altmod    : values of model altitude from 0 to 19.2km (CFMIP) or 40.5km  !
!               (LMDZ)                                                       !
!               dimension : CFMIP=41, LMDZ=19, NASA=41                       !
!               model, has dimension 41 or 19                                !
!   prestop   : same var as altmod in pressure mod                           !
!   srmod     : values of diagSR, dimension always equal to 15               !
!               values : -1,0,0.01,1.2,2,3,5,10,20,30,40,50,60,80,100        !
!   pr2moy    : observed atb averaged on altmod lvl with the routine         !
!               vertical_mean,following the satellite,                       !
!               dim=(altmax,nprof)                                           !
!   molmoy    : mol3 averaged on altmod lvl with the routine vertical_mean,  !
!               following the satellite, dim=(altmax,nprof)                  !
!   srmoy     : scattering ratio calculated with pr2moy and molmoy,          !
!               dim=(altmax,nprof)                                           !
!   indice    : number of atb iteration in an altitude box                   !
!   indicem   : number of mol iteration in an altitude box                   !
!   indiceh   : number of time iteration in a lat/lon/alt box with CHIM/WRF  !
!               model                                                        !
!   mheure    : time averaged on lon/lat/alt grid with CHIM/WRF model        !
!                                                                            !
!----------------------------------------------------------------------------!
      real*4,dimension(:),allocatable  ::  latmod, lonmod, srmod, prestop, depolmod, pr2mod, srdepmod,atbrmod, tempmod!, crmod
      real*4,dimension(:),allocatable  ::  latmid, lonmid      
      real,dimension(:),allocatable  ::  altmod, altmid,tempmid
      real,dimension(:,:),allocatable  :: altmod_bound, tempmod_bound
      real,dimension(:,:),allocatable  ::  pr2moy,  molmoy, srmoy,depolmoy, pr2moy2,crmoy,perpmoy,parmoy, tempmoy
      real,dimension(:,:),allocatable  ::  indice, indicem, indiceh, indice2,indicep,indicep2,indicetemp 
      real*4,dimension(:,:),allocatable  ::  mheure
      real*4,dimension(:,:),allocatable  :: SRwrf,CRwrf,DEPOLwrf            

!----------------------------------------------------------------------------!
! LMDZ output variables                                                      !
!----------------------------------------------------------------------------!
!               instantaneous fraction : dim=(altmax,nprof)                  !
!****************************************************************************!
!   uncertfraction   : uncertain fraction flag 0/1 for each profil and each  !
!                      altitude, calculated with routine fraction_subgrid2   !
!   satfraction      : fully attenuated flag 0/1                             !
!   cloudfraction    : cloudy fraction flag 0/1                              !
!   clearfraction    : clear fraction flag 0/1                               !
!   nanfraction      : NaN fraction flag 0/1                                 !
!   sefraction       : Surface_Elevation fraction flag 0/1                   !
!   fractot          : all 6 fraction unified in a single var, by flag 1,2,3,!
!                      4,5,6 used in the routine fraction_subgrid3 :         !
!                      1=sat, 2=clear, 3=uncert, 4=nan, 5=SE, 6=cloud        !
!                                                                            !
!----------------------------------------------------------------------------!
!               daily fraction : dim=(lonmax,latmax,altmax,day)              !

!****************************************************************************!
!   uncertfractday   : uncertfraction averaged on LMDZ/CFMIP/NASA grid day by!
!                      day.                                                  !
!   satfractday      : satfraction averaged on LMDZ/CFMIP/NASA grid          !
!   cloudfractday    : cloudfraction averaged on LMDZ/CFMIP/NASA grid        !
!   clearfractday    : clearfraction averaged on LMDZ/CFMIP/NASA grid        !
!   nanfractday      : nanfraction averaged on LMDZ/CFMIP/NASA grid          !
!   sefractday       : sefraction averaged on LMDZ/CFMIP/NASA grid           !
!   indday           : number of sat/cloud/clear/uncert values on            !
!                      lon/lat/alt/day boxes                                 !
!   inddaytot        : number of sat/cloud/clear/uncert/nan/se values on     !
!                      lon/lat/alt/day boxes                                 !
!                                                                            !
!----------------------------------------------------------------------------!
!              monthly fraction : dim=(lonmax,latmax,altmax)                 !
!****************************************************************************!
!   monthuncertfract : uncertfractday averaged on all days of the month      !
!   monthsatfract    : satfractday averaged on all days of the month         !
!   monthcloudfract  : cloudfractday on all days of the month                !
!   monthclearfract  : clearfractday on all days of the month                !
!   monthnanfract    : nanfractday averaged on all days of the month         !
!   monthsefract     : sefractday averaged on all days of the month          !
!   indpermonth      : number of sat/cloud/clear/uncert fractday on          !
!                      lon/lat/alt boxes                                     !
!   indpermonthtot   : number of sat/cloud/clear/uncert/nan/se fractday on   !
!                      lon/lat/alt boxes                                     !
!                                                                            !
!----------------------------------------------------------------------------!
!                    instantaneous isccp : dim=(nprof)                       !
!****************************************************************************!
!   isccplow        : isccp low cloud flag 0/1 for each profil               !
!   isccpmid        : isccp mid cloud flag 0/1 for each profil               !
!   isccphigh       : isccp high cloud flag 0/1 for each profil              !
!   colcloud        : isccp column cloud flag 0/1 for each profil            !
!   colclear        : isccp column clear flag 0/1 for each profil            !
!                                                                            !
!----------------------------------------------------------------------------!
!                daily isccp : dim=(lonmax,latmax,daymax)                    !

!****************************************************************************!
!   isccplowday     : isccplow averaged on LMDZ/CFMIP/NASA grid              !
!   isccpmidday     : isccpmid averaged on LMDZ/CFMIP/NASA grid              !
!   isccphighday    : isccphigh averaged on LMDZ/CFMIP/NASA grid             !
!   colcloudday     : colcloud averaged on LMDZ/CFMIP/NASA grid              !
!   colclearday     : colclear averaged on LMDZ/CFMIP/NASA grid              !
!   isccpindday     : number of isccp low/mid/high on a lon/lat/day box      !
!                                                                            !
!----------------------------------------------------------------------------!
!                      monthly isccp : dim=(lonmax,latmax)                   !
!****************************************************************************!
!   monthisccplow   : isccplowday averaged on all days of the month          !
!   monthisccpmid   : isccpmidday averaged on all days of the month          !
!   monthisccphigh  : isccphighday averaged on all days of the month         !
!   monthcolcloud   : colcloudday averaged on all days of the month          !
!   monthcolclear   : colclearday averaged on all days of the month          !
!   isccpdaypermonth: number of monthisccp low/mid/high on a lon/lat box     !
!                                                                            !
!----------------------------------------------------------------------------!
!                            diagSR variables                                !
!****************************************************************************!
!   diagSR          : number of occurence of a SR value in one of the 15     !
!                     different interval defined with srmod                  !
!                     dim=(lonmax,latmax,altmax,daymax,diagmax)              !
!   monthdiagSR     : diagSR averaged on all days of the month               !
!                     dim=(lonmax,latmax,altmax,diagmax)                     !
!   sumdiag         : sum of diagSR at an altitude                           !
!                                                                            !
!----------------------------------------------------------------------------!
      real,dimension(:,:),allocatable :: uncertfraction,satfraction,       &
                                           cloudfraction,clearfraction,      &
                                           nanfraction,sefraction,           &
                                           rejfraction,fractot,cloudfraction2

      real,dimension(:,:,:,:),allocatable :: cloudfractday, clearfractday, &
                                          
                                               uncertfractday, &
                                               indday

      real,dimension(:,:,:),allocatable :: monthcloudfract,monthclearfract,&
                                             monthuncertfract, &
                                             indpermonth,indphasepermonth

      integer,dimension(:,:,:),allocatable  ::  indnan

      integer,dimension(:),allocatable :: isccplow, isccpmid, isccphigh,      &
                                         colcloud
      integer,dimension(:,:),allocatable :: isccpliq, isccpice
      integer,dimension(:,:,:),allocatable :: isccpun

 
       integer  ::  indbase
       real,dimension(:,:,:),allocatable  ::  heightday2
       real,dimension(:),allocatable  ::  height2
       real,dimension(:,:),allocatable  ::  monthheight2

       real,dimension(:,:,:),allocatable  ::  heightday,indheight
       real,dimension(:,:),allocatable  ::  monthheight,indmonthheight
       real,dimension(:),allocatable  ::  heightmod
       integer,dimension(:,:,:),allocatable  ::  hlow,hmid,hhigh,hcol
       real,dimension(:,:,:),allocatable  ::  hheight
       real,dimension(:),allocatable  ::  height
       integer  ::  iheight, ihist, ihisttemp, icat
       integer,parameter  ::  heightmax=41, catmax=5
       integer,parameter  ::  histmax=11, histmax2=41, histtempmax=29,histtempmax2=28
       integer,dimension(histmax2)  ::  histmod2
       real,dimension(histmax)  ::  histmod
 
       real,dimension(histtempmax)  ::  histtempmod
       real,dimension(histtempmax2)  ::  histtempmod2
 
       integer,dimension(:,:,:),allocatable  :: hcoltemp,hlowtemp,hmidtemp,hhightemp
       real,dimension(:,:,:),allocatable  :: coltemp,lowtemp,midtemp,hightemp
       real,dimension(:,:,:),allocatable  :: indcoltemp,indlowtemp,indmidtemp,indhightemp
 
       real,dimension(:,:),allocatable  ::  monthlowtemp,monthmidtemp,monthhightemp,monthcoltemp
       real,dimension(:,:),allocatable  ::  indmonthlowtemp,indmonthmidtemp,indmonthhightemp,indmonthcoltemp


      real,dimension(:,:,:),allocatable :: isccplowday, isccpmidday,       &
                                             isccphighday, isccpindday,      &
                                             isccpinddaylow, isccpinddaymid, &
                                             colcloudday





      real*8  ::  colclearres

      real  ::  isccptemp

      real,dimension(:,:,:,:),allocatable :: isccpiceday, isccpliqday
      real,dimension(:,:,:,:,:),allocatable :: isccpunday
      real,dimension(:,:,:,:),allocatable :: isccpphaseday
      real,dimension(:,:,:),allocatable :: isccpdustday

      real,dimension(:,:,:),allocatable :: monthisccpice, monthisccpliq
      real,dimension(:,:,:,:),allocatable :: monthisccpun
      real,dimension(:,:,:),allocatable :: monthisccpphase,indmonthphase
      real,dimension(:,:,:),allocatable :: indmonthphase2
      integer,dimension(:,:,:),allocatable  ::  indmonthphase3D

      real,dimension(:,:),allocatable :: monthisccplow, monthisccpmid,     &
                                           monthisccphigh, isccpdaypermonth, &
                                           isccpdaypermonthlow,              &
                                           isccpdaypermonthmid,              &
                                           monthcolcloud, monthcolclear

      real,dimension(:,:,:,:),allocatable :: diagSR
      real,dimension(:,:,:,:,:),allocatable :: diagSRpha

      real,dimension(:,:,:,:,:),allocatable :: diagPHA


      integer*4  ::  sumdiag


!----------------------------------------------------------------------------!
! META variables                                                             !
!----------------------------------------------------------------------------!
!   altl   : altitude of lidar lvl in kilometer, dim=(altitude)              !
!   altm   : altitude of meteo lvl in kilometer, dim=(altitude2)             !
!                                                                            !
!----------------------------------------------------------------------------!
      real*4,dimension(:),allocatable  ::  altl, altm  
      
      integer  ::  seuilsnrlow, seuilsnrhigh
      integer  ::  nanprof, nansat, nanmid, nanlow
      integer  ::  icewaterres
      real,parameter  ::  alpha=0.0028, beta=0.0458, alpha1=3., beta1=0.0576


      real,parameter  ::  alpha50=9.0322e+03, beta50=-2.1358e+03, gamma50=173.3963, delta50=-3.9514, epsilon50=0.2559, zeta50=-9.4776e-04

      real,parameter  ::  Ahoi=0.1667, Bhoi=-0.01
      real  ::  depoltmp, perptmp1, perptmp2
      real,dimension(:,:),allocatable :: icecloud, watercloud,phasecloud
      real,dimension(:,:,:),allocatable :: uncloud
      real,dimension(:,:,:,:),allocatable :: icecloudfractday, watercloudfractday,indphaseday
      real,dimension(:,:,:,:),allocatable :: phasefractday,inddayphase
      real,dimension(:,:,:,:,:),allocatable :: uncloudfractday
      real,dimension(:,:,:,:),allocatable :: indphasefractday
      real,dimension(:,:,:,:,:),allocatable :: indphaseunday

      real,dimension(:,:,:),allocatable :: monthicecloud, monthwatercloud
      real,dimension(:,:,:),allocatable :: indphasemonth
      real,dimension(:,:,:,:),allocatable :: monthuncloud
      real,dimension(:,:,:),allocatable :: monthphasecloud
      real,dimension(:,:,:),allocatable  ::  monthcftemp,monthcftempice
      real,dimension(:,:,:),allocatable  ::  monthcftempliq
      real,dimension(:,:,:),allocatable  :: monthcftempphase,indmonthphasetemp
      real,dimension(:,:,:),allocatable  :: indcftemppermonth
      real,dimension(:,:,:,:),allocatable  :: indcftempphase,indcftemp
      real,dimension(:,:,:,:),allocatable  :: cftempday,cftempphaseday
      real,dimension(:,:,:,:),allocatable  :: cftempiceday,cftempliqday
      real,dimension(:,:),allocatable  ::  cftemp,cftempliq,cftempice
      real  ::  cfsumtemp

! variable pour fichier phase occurrences
real*4,dimension(:,:,:),allocatable  :: tot_ind,cloud_ind,ice_ind,water_ind
real*4,dimension(:,:,:,:),allocatable  :: un_ind



      metal='Lidar_Data_Altitudes'   ! name of meta var     
      metam='Met_Data_Altitudes' 



!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!
 

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                  SELECTION OF INPUT/OUPUT PARAMETER                      *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!------------------- Open output grid file -----------------------!
!                                                                 !
!   The format you have to enter is :                             !
!   "Map3D330m_",period,day/night,grid,version                    !
!   ex : "Map3D330m_200701_night_CFMIP_2.0"                       !
!                                                                 !
!-----------------------------------------------------------------!
  do
     write (*,'(a)',advance='no') 'Period identifier ? '
     read *, file5
     if (err==0) exit 
  enddo

!------------------- Open input data file ------------------------!
!                                                                 !
!   You have to enter the path and the name of the file  which    !
!   lists all the hdf CALIPSO                                     !
!   ex : "./Liste/200701_night"                       !    
!                                                                 !
!-----------------------------------------------------------------!
  do
     write (*,'(a)',advance='no') 'input file name ? '
     read *, file3
     open(unit=1,file=file3,iostat=err,status='OLD')
     if (err==0) exit
     print *,'--- input file not found' 
  enddo

!--------------------- Select the model --------------------------!
!                                                                 !
!   model = "lmdz" "chim" or "wrf" to select the LMDZ CHIMERE or  !
!           WRF model                                             !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'model (lmdz, chim or wrf) ?'
     read *, model
     if (err==0) exit 
  enddo
  
!----------------- Select the day or night -----------------------!
!                                                                 !
!   Select day or night Data version, because the process is not  !
!   the same during night and day                                 !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'night or day ? '
     read *, switch
     if (err==0) exit 
  enddo
 
!-------------------- Select the grid ----------------------------!
!                                                                 !
!    - "CFMIP" : 1° x 1° x 41levels each 480m from 0 to 19.2km    !
!             (361,181,41)                                        !
!    - "LMDZ"  : 3.75° x 2.53° x 19levels from 0 to 40.4823km     !
!             (96,72,19)                                          !
!    - "NASA"  : 5° x 5° x 41levels each 480m from 0 to 19.2km    !
!             (73,37,41)                                          !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'grid (CFMIP1, CFMIP2, LMDZ, NASA) ? '
     read *, gcm
     if (err==0) exit 
  enddo

!----------------- Select pressure or altitude -------------------!
!                                                                 !
!   Select "pressure" or "altitude" version of ouput data         !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'vertical unit (pressure, altitude) ? '
     read *, alt_pres
     if (err==0) exit 
  enddo


!----------------- Select sat or cloudy mode ---------------------!
!                                                                 !
!   Select "sat" or "cloudy" mode in order to count the first     !
!   fully attenuated point as a cloudy point in the cloudfraction !
!   if "cloudy" mode is selected, and as a fully attenuated point !
!   if "sat" mode is selected. This change appears in the routine !
!   subgrid_fraction2.                                            !
!                                                                 !
!-----------------------------------------------------------------!
  do
   write (*,'(a)',advance='no') 'sat or cloudy ? '
     read *, switch2
     if (err==0) exit 
  enddo

print *, 'input parameters entered'

!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!

!_____________________________________________________________________________
!
! Temps CPU de calcul initial.
  call cpu_time(t_cpu_0)

! Temps elapsed de reference.
  call system_clock(count=t1, count_rate=ir)
!_____________________________________________________________________________
!


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*             LOADING LON-LAT-ALT-SR-DEPOL-PR2 GRID VECTORS                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

! Read parameter file
open(5,file='./'//trim(gcm)//'.p')
   read(5,*)grid, altmax2, altmax, latmax, lonmax, toplowl, topmidl, tophighl, altfile, latfile, lonfile
 close(5)

print *, 'Grid parameter file read' 
allocate(heightmod(heightmax))
allocate(latmod(latmax),lonmod(lonmax),prestop(altmax),altmod(altmax2),srmod(diagmax),pr2mod(pr2max),atbrmod(permax), tempmod(tempmax),              &
         altmid(altmax),latmid(latmax-1),lonmid(lonmax-1), depolmod(depolmax), srdepmod(pr2max),altmod_bound(altmax,2))
allocate(tempmod_bound(tempmax-1,2),tempmid(tempmax-1))
  heightmod(:)=0;
  prestop(:)=0;latmod(:)=0;lonmod(:)=0;
  altmod(:)=0;srmod(:)=0; depolmod(:)=0; pr2mod(:)=0;atbrmod(:)=0;tempmod(:)=0;
  altmid(:)=0; latmid(:)=0; lonmid(:)=0;
  srdepmod(:)=0;

altmod_bound(:,:)=0; tempmod_bound(:,:)=0;tempmid(:)=0;
  
! loading the grid of the diagSR boxes value
open(17,file='./grilles_lmdz/srmod10')
!print *, 'open the file'

  do idiag=1,diagmax
     read(17,*)srmod(idiag)
  enddo
 close(17)

! loading the grid of the DepolSR boxes value
open(6,file='./grilles_lmdz/depolmod')

  do idep=1,depolmax
     read(6,*)depolmod(idep)
  enddo
 close(6)

! loading the grid of the pr2 boxes value
open(18,file='./grilles_lmdz/atbmod301')
  do ipr2=1,pr2max
     read(18,*)pr2mod(ipr2)
  enddo
 close(18)

open(23,file='./grilles_lmdz/atbrmod241')
  do ipr2=1,permax
     read(23,*)atbrmod(ipr2)
  enddo
 close(23)

! loading the grid of the DepolSR boxes value
open(20,file='./grilles_lmdz/tempmod39')

  do itemp=1,tempmax
     read(20,*)tempmod(itemp)
   enddo
 close(20)

        do itemp=1,tempmax-1
            tempmid(itemp) = (tempmod(itemp)+tempmod(itemp+1))/2
            tempmod_bound(itemp,1)=tempmod(itemp);
            tempmod_bound(itemp,2)=tempmod(itemp+1);
         enddo

! Computing the Height grid
 do iheight=1,heightmax-1
    heightmod(iheight+1)=heightmod(iheight)+0.5
 enddo



! loading the level grid (altitude or pressure)
    if(alt_pres=='altitude')then

open(15,file='./grilles_lmdz/'//altfile)
         do iz=1,altmax2
            read(15,*)altmod(iz)
         enddo
         do iz=1,altmax
            altmid(iz) = (altmod(iz)+altmod(iz+1))/2
            altmod_bound(iz,1)=altmod(iz);
            altmod_bound(iz,2)=altmod(iz+1);
         enddo

         do iz=1,altmax     
            if(altmod(iz+1).GE.8.64)then
               seuilsnrhigh=iz
               exit
            endif
         enddo

         do iz=altmax,1,-1     
            if(altmod(iz+1).LE.2.4)then
                seuilsnrlow=iz
                exit
            endif
         enddo

    elseif(alt_pres=='pressure')then
         if(trim(gcm).eq.'LMDZ')then
            open(15,file='./grilles_lmdz/pression_lmdz2.txt')
         else if(trim(gcm).eq.'CFMIP')then
            open(15,file='./grilles_lmdz/pres_cfmip')
         elseif(trim(gcm).eq.'NASA')then
            open(15,file='./grilles_lmdz/pres_cfmip')   

         endif
            do iz=1,altmax
            read(15,*)prestop(iz)           ! lmdz milieu de la couche 
            enddo       
   endif
 close(15)


open(10,file='./grilles_lmdz/'//lonfile)
      do ilon=1,lonmax
         read(10,*)lonmod(ilon)           !lmdz
      enddo
      do ilon=1,lonmax-1
         lonmid(ilon)=(lonmod(ilon)+lonmod(ilon+1))/2
      enddo
 close(10)

open(21,file='./grilles_lmdz/'//latfile)
      do ilat=1,latmax
         read(21,*)latmod(ilat)           ! lmdz
      enddo
      do ilat=1,latmax-1
         latmid(ilat)=(latmod(ilat)+latmod(ilat+1))/2
      enddo
 close(21)
  print *, 'latmod lonmod altmod srmod depolmod pr2mod ok'

!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!



!____________________________________________________________________________
!
! Initialization of box and numfich

      box=0
      numfich=0

! Reading the list of hdf files

 888  read(1,100,end=999)file2

      print *, 'Processing with ',trim(file2)
      numfich=numfich+1
      print *, 'lecture du fichier numero ',numfich

! initialization of var
   it=0; comptpf=0; date=0; year=0; month=0; day=0; 
filetmp2=trim(file2)



!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*             BEGINING OF THE READING OF THE SDS/META VAR                  *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

print *, 'Read the Calipso file data ...'


!****************************** READING SDS VAR ****************************!

 call nprof(filetmp2,20,it)       ! find the number of profil it

! empty file checking
if(it.lt.500)then
goto 887
endif


!indtot=indtot+(it*ilid)            ! indice total of data used

! Allocation of interpolated variables
allocate(lat(it),lon(it),SE(it),temps(it),temps2(it), &
          mol2(altitude,it),   & 
         mol3(altitude, it), temp2(altitude, it), stat = OK_buffer)


if(alt_pres=='pressure')then
   allocate(pres2(altitude, it))
   pres2(:,:)=0
endif

! check the allocation 
 if (OK_buffer/=0) print *,'--- buffer allocation error '   

! Initialization of interpolated variables
   temps2(:)=0; temps(:)=0; mol2(:,:)=0; mol3(:,:)=0; temp2(:,:)=0;
   lat(:)=0;lon(:)=0;SE(:)=0;




!---------------- Select the type of hdf file ( prov/launch )----------------!
!                                                                            !
! Depending on the version of the hdf files (prov or launch), the time       !
! variable is different. That's why, the process isn't the same in the 2cases!
!                                                                            !
!----------------------------------------------------------------------------!

if(file2(67:70)/="Laun")then

   !     Retrieve data for time variable.
   sds_varname='Profile_UTC_Time'
   call sdsread8(vartmp8,filetmp2,sds_varname,ret1)
   
   if(ret1.eq.-1)then
      deallocate(vartmp8)
      deallocate(lat,lon,SE,temps,temps2,mol2, mol3, mol4,indtot, indret, stat = OK_buffer)
      deallocate(indretlow,indretmid,indrethigh)
      command3='echo '//trim(file2)//' >> ./cal_corrompu'
      call system(trim(command3))
      goto 887
   endif
   	
   temps(:)=vartmp8(1,:);
   deallocate(vartmp8)


!------------------ Calculation of the date : type Prov ---------------------!
!                                                                            !
!   temps(nprof) is in International Atomic Time in UTC : yymmdd.ffffffff    !
!                                                                            !
!   yy  =  Last two digits of year where 07 represents 2007                  !
!   mm  =  Month in two-character subfield with values 01-12                 !
!   dd  =  Day of month in two-character subfield with values 01-28, -29,..  !
!   "."  =  Period as a separator                                            !
!   ffffffff  =  Fractional part of day                                      !
!                                                                            !
!----------------------------------------------------------------------------!

            date=int(temps(1))           ! date of the read file
            year=2000+int(date/10000)                      ! year
            month=int((date-int(date/10000)*10000)/100)    ! month
            day=date-int(date/10000)*10000-month*100       ! day
            jour=day
            print *, 'Processing for the ',day,'/',month,'/',year
            print *, 'File type = Prov or ValStage'

else 

!------------ Calculation of the date : type Launch (old data) --------------!
!                                                                            !
! If the file is a launch one then, the UTC time is calculated from the TAI  !
! time in seconds.                                                           !
! This calculation is necessary for the June July and August months period.  !
!                                                                            !
!----------------------------------------------------------------------------!

   sds_varname='Profile_Time'
   call sdsread8(vartmp8,filetmp2,sds_varname,ret1)
   temps(:)=vartmp8(1,:);
   deallocate(vartmp8)


   date=int(temps(1)) 
   year=2006

   if(date.lt.425865606)then   ! first value of the time in June
      month=6
      day=int((date-423273606)/86400)+1  
   elseif(date.lt.428544006)then   ! first value of time in July
         month=7
         day=int((date-425865606)/86400)+1
   else
         month=8
         day=int((date-428544006)/86400)+1   ! first value of time in August
   endif

   print *, 'Processing for the ',day,'/',month,'/',year
   print *, 'Fichier Launch'

endif ! end selection of file type


!---------------- Time in day for the output netcdf files -------------------!
!                                                                            !
! Select the number of days since 2000-01-01 00:00:00 to the mid of the run  !
! period.                                                                    !
!                                                                            !
! ex : for monthly period, date is set to the 15th of the month              !
!      for a trimonthly period date is set to the 15th of the second month   !
!                                                                            !
!----------------------------------------------------------------------------!

if(numfich.eq.1)then 
     print *, file3(25:30)
     date1=((year*100+month)*100+15)*100
     date2=2000010100
     print *, date1
     print *, date2

     resh=interdat(date2,date1) ! calculation of time between the 2 period 
     print *, resh              ! in hour
     resd=resh/24               ! converted in days
     print *, resd
endif


!     Retrieving data for atb variable.
sds_varname='Total_Attenuated_Backscatter_532'
 call sdsread(atb,filetmp2,sds_varname)

!     Retrieving data for atb variable.
sds_varname='Attenuated_Backscatter_1064'
 call sdsread(atb2,filetmp2,sds_varname)

!     Retrieving data for atb variable.
sds_varname='Perpendicular_Attenuated_Backscatter_532'
 call sdsread(perp,filetmp2,sds_varname)

!     Retrieving data for lat variable.
sds_varname='Latitude'
 call sdsread(vartmp,filetmp2,sds_varname)
lat(:)=vartmp(1,:);
deallocate(vartmp)

!     Retrieving data for lon variable.
sds_varname='Longitude'
 call sdsread(vartmp,filetmp2,sds_varname)
lon(:)=vartmp(1,:);
deallocate(vartmp)

!     Retrieving data for mol variable.
sds_varname='Molecular_Number_Density'
 call sdsread(mol,filetmp2,sds_varname) 

!     Retrieve data for temp variable.
sds_varname='Temperature'
 call sdsread(temp,filetmp2,sds_varname)

!     Retrieving data for surf_elevation !!
sds_varname='Surface_Elevation'
 call sdsread(vartmp,filetmp2,sds_varname)
SE(:)=vartmp(1,:);
deallocate(vartmp)


!***************************** READING META VAR *****************************!

!     Retrieving data for altitude variable.
 call metaread(altl,metal,filetmp2)
!     Retrieving data for altitude variable.
 call metaread(altm,metam,filetmp2)


print *, 'HDF Calipso File read'
print *, 'Input variables read'
print *, ''
!_____________________________________________________________________________
!



print *, 'Interpolation of data & molecular calculation start ...'

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*           BEGINNING OF CALCULATION ON THE INPUT VARIABLES                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!



!************** CONVERTION OF THE UTC FRACTION IN UT HOUR *******************!



!******************************* WRF  MODE *********************************!
if(model.eq.'wrf')then

do i=1,it          !!!! loop on each profil   
   if((file2(67:70)=="Prov").or.(file2(67:70)=="AL_L"))then
   ! Calculation of hour UT with UTC time fraction
      temps2(i)=(temps(i)-date)*24  ! temps = fraction of hour in prov file type
   else 
   ! Calculation of hour UT with TAI time
      if(file2(67:70)=="Laun")then
         if(date.lt.425865606)then
            temps2(i)=(((temps(i)-423273606)/86400+1)-day)*24
         elseif(date.lt.428544006)then
            temps2(i)=(((temps(i)-425865606)/86400+1)-day)*24    
         else
            temps2(i)=(((temps(i)-428544006)/86400+1)-day)*24  
         endif
      endif
   endif
enddo

endif
!_____________________________________________________________________________
!


!****************************** LAUNCH FILE *********************************!   
if(file2(67:70)=="Laun")then
   do i=1,it          !!!! loop on each profil
      jour=day
      date=60000+day+100*month
   enddo
endif
!_____________________________________________________________________________
!


!  INTERPOLATION OF PRESSURE TEMPERATURE AND MOLECULAR FROM 33LVL TO 583LVL  !
if(alt_pres=='pressure')then
 do i=1,it          !!!! loop on each profil
   call interp(pres,pres2,altm,altl,i,it)
   call interp(mol,mol2,altm,altl,i,it)     ! atb molecular without dimension

 enddo        !!!! end of profil loop
else
 do i=1,it          !!!! loop on each profil
   temps2(i)=(temps(i)-date)*24  ! temps = fraction of hour in prov file type
   call interp(mol,mol2,altm,altl,i,it)     ! atb molecular without dimension
   call interp(temp,temp2,altm,altl,i,it)
   call atb_temp_interp(temp2,altl,i,it,SeuilTemp1km,SE)         ! molecular extrapolation
   call SE_alt_mol(SE,altl,temp2,altitude,it,i)! add the SE to mol3

enddo        !!!! end of profil loop
endif


deallocate(temps)


!************************ CALCULATION OF ATB-MOLECULAR **********************!
 

!---------------------------- day or night mode -----------------------------!
!                                                                            !
!  The calculation of normalized ratio is calculated at different altitude   !
!  levels whether it is day or night time.                                   !
!  This difference appears in the "atb_mol" routine:                         !
!     - between 20km and 25km during day time                                !
!     - between 22km and 25km during night time                              !
!                                                                            !
!----------------------------------------------------------------------------!

!-------------------------------- PSC filter --------------------------------!
!                                                                            !
!  During the polar winter, the cloud could appear up to 20km altitude at    !
!  latitude higher than 60° south. That's why from June to October, the      !
!  altitude calculation of the normalized ratio must be higher at these      !
!  latitudes.                                                                !
!  The calculation altitude is ranged :                                      !
!     - between 28.5km and 35km during day time                              !
!     - between 28.5km and 33km during night time                            !
!                                                                            !
!----------------------------------------------------------------------------!
if(switch.eq.'day')then   !!! DAY mode  !!!

 do i=1,it          !!!! loop on each profil   
    ! PSC filter during the polar winter
   if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then 
      ! normalized ratio betwen 28.5 & 35km & molecr calculation in km-1 sr-1  
      call atb_mol(atb,mol2,mol3,i,it,17,42) ! Antarctic PSCs season 
   elseif((lat(i).ge.60).and.((month.eq.12).or.(month.eq.1).or.(month.eq.2)))then
      call atb_mol(atb,mol2,mol3,i,it,17,42) ! Arctic PSCs season   Pitts et al 2011   
   else
      ! normalized ratio betwen 20 & 25km & molecular calculation in km-1 sr-1
      call atb_mol(atb,mol2,mol3,i,it,62,92) 
   endif



      call atb_mol_interp(mol3,altl,i,it,SeuilMol1km,SE)         ! molecular extrapolation
      call SE_alt_mol(SE,altl,mol3,altitude,it,i)! add the SE to mol3
      call SE_alt_atb(SE,altl,atb,altitude,it,i) ! add the SE to atb

 enddo         !!!! end of profil loop

elseif(switch.eq.'night')then  !!! NIGHT mode !!!

 do i=1,it          !!!! loop on each profil  
  ! PSC filter during the polar winter
  if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
     ! normalized ratio betwen 28.5 & 33km & molecr calculation in km-1 sr-1  
     call atb_mol(atb,mol2,mol3,i,it,24,42)  
  elseif((lat(i).ge.60).and.((month.eq.12).or.(month.eq.1).or.(month.eq.2)))then
      call atb_mol(atb,mol2,mol3,i,it,24,42) ! Arctic PSCs season   Pitts et al 2011
  else
     ! normalized ratio betwen 22& 25km & molecular calculation in km-1 sr-1
     call atb_mol(atb,mol2,mol3,i,it,62,78)  
  endif


      call atb_mol_interp(mol3,altl,i,it,SeuilMol1km,SE)        ! molecular extrapolation
      call SE_alt_mol(SE,altl,mol3,altitude,it,i) ! add the SE to mol3
      call SE_alt_atb(SE,altl,atb,altitude,it,i)  ! add the SE to atb
 enddo         !!!! end of profil loop


endif  ! end of day/night loop

print *, 'Interpolation of data & molecular calculation done'
!_____________________________________________________________________________
!




!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*            SELECTION OF THE MODEL GRID CFMIP / LMDZ / NASA               *!
!*               AND ALLOCATION OF THE OUTPUTS VARIABLES                    *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

    
! noms output : base + period + day/night + grille + sat/cloudy + version
! e.g. 
! SR_histo330m_201005_night_CFMIP2_sat_2.68.nc
! base = SR_histo330m_
! period = 200701 [on va recycler file5 pour ca]
! day/night  [variable switch]
! grille = LMDZ, CFMIP1, CFMIP2, NASA [variable gcm]
! sat/cloudy [variable switch2]
! version = 2.68 [variable version]:

    
! Output Map file name
file6='MapLowMidHigh_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
! Output Map file name
file66='MapHigh_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
! Output diagSR file name
file7='SR_histo_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
! Output diagSR file name
file12='SR_histo_Phase_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
! Output depolSR file name
file10='3D_CloudFraction_Phase_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
file13='3D_CloudFraction_Temp_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
! Output Map file name
file11='MapLowMidHigh_Phase_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)
! Output Map file name
! file12='Phase_histo'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)

print *, file6
print *, 'altmax=',altmax


! Allocation / initialization of variable 
print *, 'Allocation / Initialization of variables...'

  allocate(indicep(altmax,it),indicep2(altmax,it),pr2moy2(altmax,it),indice2(altmax,it),crmoy(altmax,it))
  allocate(pr2moy(altmax,it),molmoy(altmax,it),srmoy(altmax,it),depolmoy(altmax,it),perpmoy(altmax,it),parmoy(altmax,it),tempmoy(altmax,it))
  allocate(indice(altmax,it),indicem(altmax,it),indicetemp(altmax,it)) 


if(numfich.eq.1)then ! Allocation of monthly variables


   allocate(indday(latmax-1,lonmax-1,altmax,daymax))

   allocate(cloudfractday(latmax-1,lonmax-1,altmax,daymax),                      &
            clearfractday(latmax-1,lonmax-1,altmax,daymax))

   allocate(uncertfractday(latmax-1,lonmax-1,altmax,daymax))

   allocate(icecloudfractday(latmax-1,lonmax-1,altmax,daymax),               &
            phasefractday(latmax-1,lonmax-1,altmax,daymax),                      &
            watercloudfractday(latmax-1,lonmax-1,altmax,daymax))
   allocate(uncloudfractday(latmax-1,lonmax-1,altmax,daymax,catmax))

   allocate(indphaseday(latmax-1,lonmax-1,altmax,daymax))
   allocate(indphaseunday(latmax-1,lonmax-1,altmax,daymax,catmax))

   allocate(isccplowday(latmax-1,lonmax-1,daymax),                               &
            isccpmidday(latmax-1,lonmax-1,daymax),                               &
            isccphighday(latmax-1,lonmax-1,daymax))

   allocate(heightday(latmax-1,lonmax-1,daymax),                                &
            indheight(latmax-1,lonmax-1,daymax))
   allocate(heightday2(latmax-1,lonmax-1,daymax))

   allocate(colcloudday(latmax-1,lonmax-1,daymax),                               &

            isccpinddaylow(latmax-1,lonmax-1,daymax),                            &
            isccpinddaymid(latmax-1,lonmax-1,daymax),                            &
            isccpindday(latmax-1,lonmax-1,daymax))


allocate(cftempday(latmax-1,lonmax-1,tempmax-1,daymax), &
         cftempiceday(latmax-1,lonmax-1,tempmax-1,daymax), &
         cftempliqday(latmax-1,lonmax-1,tempmax-1,daymax), &
         indcftemp(latmax-1,lonmax-1,tempmax-1,daymax), &
         indcftempphase(latmax-1,lonmax-1,tempmax-1,daymax))

   allocate(lowtemp(latmax-1,lonmax-1,daymax),                         &
            midtemp(latmax-1,lonmax-1,daymax),                          &
            hightemp(latmax-1,lonmax-1,daymax),                          &
            coltemp(latmax-1,lonmax-1,daymax))
   allocate(indlowtemp(latmax-1,lonmax-1,daymax),                         &
            indmidtemp(latmax-1,lonmax-1,daymax),                          &
            indhightemp(latmax-1,lonmax-1,daymax),                          &
            indcoltemp(latmax-1,lonmax-1,daymax))



   allocate(diagSR(lonmax-1,latmax-1,altmax,diagmax-1))
   allocate(diagSRpha(lonmax-1,latmax-1,altmax,diagmax-8,3))

   allocate(diagPHA(lonmax-1,latmax-1,altmax,tempmax-1,3))

   allocate(indnan(latmax-1,lonmax-1,altmax))


allocate(indtotmean(lonmax-1,latmax-1),indtot(lonmax-1,latmax-1))
   allocate(isccpliqday(latmax-1,lonmax-1,daymax,4),isccpiceday(latmax-1,lonmax-1,daymax,4))
   allocate(isccpunday(latmax-1,lonmax-1,daymax,4,catmax),                               &
            isccpphaseday(latmax-1,lonmax-1,daymax,4))

endif

if(alt_pres=='pressure')then
   allocate(SEp(it))
   SEp(:)=0
endif

  pr2moy(:,:)=0; molmoy(:,:)=0;tempmoy(:,:)=0;  
  srmoy(:,:)=0; depolmoy(:,:)=0;
  indice(:,:)=0; indicem(:,:)=0;
  pr2moy2(:,:)=0;indice2(:,:)=0;crmoy(:,:)=0;
  indicep2(:,:)=0; indicep(:,:)=0; indicetemp(:,:)=0;

if(numfich.eq.1)then 
   indday(:,:,:,:)=0;
  cloudfractday(:,:,:,:)=0;clearfractday(:,:,:,:)=0;
  uncertfractday(:,:,:,:)=0;
  icecloudfractday(:,:,:,:)=0; watercloudfractday(:,:,:,:)=0
  uncloudfractday(:,:,:,:,:)=0;
  phasefractday(:,:,:,:)=0; 

indcftempphase(:,:,:,:)=0
indcftemp(:,:,:,:)=0
 cftempday(:,:,:,:)=0
 cftempiceday(:,:,:,:)=0
 cftempliqday(:,:,:,:)=0

  indphaseday(:,:,:,:)=0; indphaseunday(:,:,:,:,:)=0; 
  heightday(:,:,:)=0; indheight(:,:,:)=0;
  heightday2(:,:,:)=0;
  isccplowday(:,:,:)=0;isccpmidday(:,:,:)=0;isccphighday(:,:,:)=0;
  colcloudday(:,:,:)=0;isccpindday(:,:,:)=0;
  isccpinddaylow(:,:,:)=0;isccpinddaymid(:,:,:)=0;
  diagSR(:,:,:,:)=0;
  diagSRpha(:,:,:,:,:)=0;

  diagPHA(:,:,:,:,:)=0;

  isccpliqday(:,:,:,:)=0;isccpiceday(:,:,:,:)=0;
  isccpunday(:,:,:,:,:)=0; isccpphaseday(:,:,:,:)=0;

lowtemp(:,:,:)=0;midtemp(:,:,:)=0;hightemp(:,:,:)=0;coltemp(:,:,:)=0;
indlowtemp(:,:,:)=0;indmidtemp(:,:,:)=0;indhightemp(:,:,:)=0;indcoltemp(:,:,:)=0;


indtotmean(:,:)=0; indtot(:,:)=0;

endif
print *, ''
!_____________________________________________________________________________
!


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*       AVERAGING OF OBSERVATIONS DATA OVER THE VERTICAL MODEL BOXES       *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


print *, 'Begin of vertical average'

!**************** VERTICAL AVERAGE OF OBSERVATIONS DATA *********************!



if(alt_pres=='altitude')then


   do i=1, it         !!! BEGIN OF IT LOOP
      do iz=altmax,1,-1 
         do ilid=1,altitude
            if ( (altl(ilid).ge.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
            then   

             call vertical_mean(temp2,tempmoy,atb,indicetemp,i,iz,ilid,it,altmax,&
                                 altitude,4)    
              call vertical_mean(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,  &
                                  altitude,1)      
              call vertical_mean(mol3,molmoy,atb,indicem,i,iz,ilid,it,altmax,&
                                  altitude,3)
                call vertical_mean(atb2,pr2moy2,atb,indice2,i,iz,ilid,it,altmax,&
                                  altitude,3)    
               call vertical_mean(perp,perpmoy,atb,indicep,i,iz,ilid,it,altmax,&
                                  altitude,3)    
               call vertical_mean(perp,parmoy,atb,indicep2,i,iz,ilid,it,altmax,&
                                  altitude,2)                   

            endif
         enddo
      enddo

      do iz=altmax,1,-1 
         call zero_detect(pr2moy,i,iz,it,altmax)
         call zero_detect(molmoy,i,iz,it,altmax)
         call zero_detect(tempmoy,i,iz,it,altmax)  
         call zero_detect(pr2moy2,i,iz,it,altmax)
         call zero_detect(parmoy,i,iz,it,altmax)
         call zero_detect(perpmoy,i,iz,it,altmax)
      enddo

          call Surf_detect2(SE,altmod,pr2moy,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,molmoy,altmax,it,i,alt_pres)  
          call Surf_detect2(SE,altmod,tempmoy,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,pr2moy2,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,parmoy,altmax,it,i,alt_pres)
          call Surf_detect2(SE,altmod,perpmoy,altmax,it,i,alt_pres)

     do iz=altmax,1,-1 
       if((trim(gcm).eq.'LMDZ').or.(trim(gcm).eq.'WRF'))then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
            call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
            call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                                 altmax)
            call SR_CR_DEPOL_mean(parmoy,perpmoy,depolmoy,indicep2,indicep,i,iz,it,   &
                                 altmax)
         else
            call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
            call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                                 altmax)
            call SR_CR_DEPOL_mean(perpmoy,parmoy,depolmoy,indicep,indicep2,i,iz,it,   &
                                 altmax)
         endif  
       
       elseif(trim(gcm).eq.'LMDZ40')then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
           call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
         else
           call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
         endif
       
       elseif((trim(gcm).eq.'CFMIP').or.(trim(gcm).eq.'CFMIP1').or.(trim(gcm)&
       .eq.'CFMIP2.5').or.(trim(gcm).eq.'CFMIP80').or.(trim(gcm).eq.'CFMIP2').or.(trim(gcm).eq.'CFMIP160').or.(trim(gcm).eq.'CFMIP320'))then


         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(perpmoy,parmoy,depolmoy,indicep2,indicep,i,iz,it,   &
                               altmax)


       elseif(trim(gcm).eq.'NASA')then 
         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
     
       endif
     enddo

   !    Detection of the Surface Elevation & set -888 for value under this 
   !    threshold

      call Surf_detect(SE,altmod,srmoy,altmax,it,i)
      call Surf_detect(SE,altmod,crmoy,altmax,it,i)
      call Surf_detect(SE,altmod,depolmoy,altmax,it,i)
      

   enddo           !!! END OF IT LOOP


elseif(alt_pres=='pressure')then
   
do i=1, it         !!! BEGIN OF IT LOOP
      do iz=altmax,1,-1
         do ilid=1,altitude
         if ( (pres2(ilid,i).gt.prestop(iz)).and.(pres2(ilid,i).lt.          &
            prestop(iz-1)) )then
              call vertical_mean(atb,pr2moy,atb,indice,i,iz,ilid,it,altmax,  &
                                  altitude,1)     
              call vertical_mean(mol3,molmoy,atb,indicem,i,iz,ilid,it,altmax,&
                                  altitude,3)    
            endif
         enddo
      enddo
      do iz=altmax,1,-1
         call zero_detect(pr2moy,i,iz,it,altmax)
         call zero_detect(molmoy,i,iz,it,altmax)
       
      enddo


         call SE_km_2_Pres2(SE,SEp,altl,pres2,prestop,pr2moy,altmax,it,i) 
         call SE_km_2_Pres2(SE,SEp,altl,pres2,prestop,molmoy,altmax,it,i) 

     do iz=altmax,1,-1 
       if((trim(gcm).eq.'LMDZ').or.(trim(gcm).eq.'WRF'))then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
            call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
      
         else
            call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
      
         endif  
       
       elseif(trim(gcm).eq.'LMDZ40')then
         if((lat(i).le.-60).and.(month.ge.6).and.(month.le.10))then
           call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it, &
                                 altmax)
         else
           call filtre_2lvl(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,      &
                            altmax,gcm)
         endif
       
       elseif((trim(gcm).eq.'CFMIP').or.(trim(gcm).eq.'CFMIP1').or.(trim(gcm)&
       .eq.'CFMIP2.5').or.(trim(gcm).eq.'CFMIP80').or.(trim(gcm).eq.'CFMIP2').or.(trim(gcm).eq.'CFMIP160').or.(trim(gcm).eq.'CFMIP320'))then


    if ((pr2moy(iz,i).eq.(-9999)).or.(molmoy(iz,i).eq.(-9999))) &
       then
       continue
     else


         do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indret(i)=indret(i)+1
                     endif
             endif
          enddo




        if (iz.lt.toplowl) then
          do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indretlow(i)=indretlow(i)+1
                     endif
             endif
          enddo
       endif

       if ((iz.ge.toplowl).and.(iz.le.topmidl))then
          do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indretmid(i)=indretmid(i)+1
                     endif
             endif
          enddo
       endif

       if(iz.gt.topmidl) then
          do ilid=1,altitude
             if ( (altl(ilid).gt.altmod(iz)).and.(altl(ilid).lt.altmod(iz+1)))&
                  then  
                     if(mol4(ilid,i).ne.-9999)then
                        indrethigh(i)=indrethigh(i)+1
                     endif
             endif
          enddo
       endif

       test=indretlow(i)+indretmid(i)+indrethigh(i)
       if(test.ne.indret(i))print *, 'error indice retrieval', indret

    endif

         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(pr2moy2,pr2moy,crmoy,indice2,indice,i,iz,it,   &
                               altmax)
         call SR_CR_DEPOL_mean(parmoy,perpmoy,depolmoy,indicep2,indicep,i,iz,it,   &
                               altmax)

       elseif(trim(gcm).eq.'NASA')then 
         call SR_CR_DEPOL_mean(pr2moy,molmoy,srmoy,indice,indicem,i,iz,it,   &
                               altmax)
        
       endif
     enddo


      call SE_km_2_Pres(SE,SEp,altl,pres2,prestop,srmoy,altmax,it,i) 
 enddo
endif            !!! END OF IT LOOP


! Mode wrf model
if(model.eq.'wrf')goto 622

print *, 'deallocate input var'


! Deallocate SDS/META variables
deallocate(atb,stat = OK_buffer)
deallocate(atb2,stat = OK_buffer)
deallocate(perp,stat = OK_buffer)
deallocate(mol,stat = OK_buffer)
!deallocate(pres,stat = OK_buffer) 
deallocate(mol2,stat = OK_buffer)
deallocate(mol3,stat = OK_buffer)
deallocate(temp,stat = OK_buffer)
deallocate(temp2,stat = OK_buffer)


!deallocate(mol4)
if(alt_pres=='pressure')deallocate(pres2,stat = OK_buffer)
deallocate(altl,stat = OK_buffer)
deallocate(altm,stat = OK_buffer)



!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                 CLOUD DIAGNOSTICS FOR ONE PROFIL                         *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

print *, 'allocate flag var'
! Allocation / initialization of instantaneous fraction variables
   allocate(cloudfraction(altmax,it),clearfraction(altmax,it),satfraction    &
            (altmax,it), uncertfraction(altmax,it),rejfraction(altmax,it))
   allocate(nanfraction(altmax,it),sefraction(altmax,it))
uncertfraction(:,:)=0; nanfraction(:,:)=0;sefraction(:,:)=0;
rejfraction(:,:)=0;



!****************************************************************************!
!******************** INSTANTANEOUS FRACTION DIAGNOSTIC *********************! 
!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                   RECORD INSTANTANEOUS SR FILES                          *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

seuilatb=2.5e-03
deltatb=0

! select the instant file type: 
!  - on = instant classic file activated
!  - off = desactivated
!  - fraction = instant fraction file with GOCCP cloud mask

! FIXME : this could be commandline argument

instant_switch='on'

SELECT CASE(instant_switch)


!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
!*                               RECORD ON                                  *!
!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
! instant classic file with SR value
CASE('on')

print *, 'instant_switch = on'

allocate(cloudfraction2(altmax,it))
cloudfraction2(:,:)=0;

write(numfichc,'(i4)')numfich
write(datec,'(I6)')date
write(yearc,'(I4)')year

! FIXME : could be better
! extracts orbit identifier, i.e.
! translates /bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.01/2007/2007_01_01/CAL_LID_L1-ValStage 1-V3-01.2007-01-01T00-22-49ZN.hdf in 2007-01-01T00-22-49ZN
command4='echo '//trim(file2)//'| cut -d/ -f8 | cut -d. -f2 > ./instant/instantname'//yearc//datec(3:6)//'_'//trim(switch)//'_'//trim(gcm)

CALL SYSTEM(trim(command4))
open(10,file='./instant/instantname'//yearc//datec(3:6)//'_'//trim(switch)//'_'//trim(gcm))
read(10,*)instantname
 close(10)

do i=1,it       !!!!! BEGIN OF IT LOOP 
! flag 0/1        

  call fraction_subgrid2_8km(seuilsnrlow,seuilsnrhigh,srmoy,pr2moy,indice,   &
                             molmoy,indicem,satfraction, &
                             cloudfraction,clearfraction,uncertfraction,     &
                             nanfraction,sefraction,rejfraction,i,altmax,it, &
                             toplowl,topmidl,switch,switch2) 
  call fraction_subgrid3_8km(seuilsnrlow,seuilsnrhigh,altmod,srmoy,pr2moy,     &
                             indice,molmoy,indicem,cloudfraction2,i,altmax,it,  &
                             switch,switch2)  

enddo      !!! END OF IT LOOP   


where(rejfraction.eq.1)
srmoy=-777.
 crmoy=-777.
depolmoy=-777.
pr2moy=-777.
molmoy=-777.
parmoy=-777.
perpmoy=-777.
endwhere

where(nanfraction.eq.1)
srmoy=-9999.
 crmoy=-9999.
depolmoy=-9999.
pr2moy=-9999.
molmoy=-9999.
parmoy=-9999.
perpmoy=-9999.

endwhere

where(sefraction.eq.1)
srmoy=-888.
 crmoy=-888.
depolmoy=-888.
pr2moy=-888.
molmoy=-888.
parmoy=-888.
perpmoy=-888.
tempmoy=-888.
endwhere

where((srmoy.ne.-9999.).and.(srmoy.ne.-777.).and.(srmoy.ne.-888.))
pr2moy=pr2moy/indice
molmoy= molmoy/indicem
parmoy= parmoy/indicep2
perpmoy= perpmoy/indicep
endwhere

where(tempmoy.ne.-888.)
tempmoy= tempmoy/indicetemp
endwhere

!!!!! RECORD INSTANT SR FILES WITH ATB ATBper ATBpar ATBmol
file4='./instant/instant_SR_CR_DR_'//trim(instantname)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(version)//'.nc'
print *, 'recording instant SR CR DR file : ', file4
 call SR_CR_DR_ATB_nc(file4,altmid,altmod_bound,resd,altmax,switch,gcm,it,lat,lon,SE,temps2,&
                  srmoy,crmoy,depolmoy,pr2moy,molmoy,perpmoy,parmoy,tempmoy,cloudfraction2)


deallocate(cloudfraction2)


!****************************************************************************!
!****************************************************************************!
!****************************************************************************!
!*                               RECORD OFF                                 *!
!****************************************************************************!
!****************************************************************************!
!****************************************************************************!

! Do not record any instant file
CASE('off')
 
print *, 'instant_switch = off'
 
do i=1,it       !!!!! BEGIN OF IT LOOP 

  call fraction_subgrid2_8km(seuilsnrlow,seuilsnrhigh,srmoy,pr2moy,indice,   &
                             molmoy,indicem,satfraction, &
                             cloudfraction,clearfraction,uncertfraction,     &
                             nanfraction,sefraction,rejfraction,i,altmax,it, &
                             toplowl,topmidl,switch,switch2) 


!*************** instant SR corrected by delta atb ****************!


enddo      !!! END OF IT LOOP   


where(rejfraction.eq.1)
srmoy=-777.
crmoy=-777.
depolmoy=-777.
pr2moy=-777.
molmoy=-777.
parmoy=-777.
perpmoy=-777.

endwhere

where(nanfraction.eq.1)
srmoy=-9999.
 crmoy=-9999.
depolmoy=-9999.
pr2moy=-9999.
molmoy=-9999.
parmoy=-9999.
perpmoy=-9999.

endwhere

where(sefraction.eq.1)
srmoy=-888.
 crmoy=-888.
depolmoy=-888.
pr2moy=-888.
molmoy=-888.
parmoy=-888.
perpmoy=-888.
tempmoy=-888.
endwhere

where((srmoy.ne.-9999.).and.(srmoy.ne.-777.).and.(srmoy.ne.-888.))
pr2moy=pr2moy/indice
molmoy= molmoy/indicem
parmoy= parmoy/indicep2
perpmoy= perpmoy/indicep
endwhere

where(tempmoy.ne.-888.)
tempmoy= tempmoy/indicetemp
endwhere

CASE DEFAULT

print *, 'instant_switch is invalid : ', instant_switch


END SELECT

! End of instant_switch select case


if(allocated(temps2)) deallocate(temps2,stat = OK_buffer)
deallocate(SE,stat = OK_buffer)



!************** INSTANTANEOUS ISCCP LOW MID HIGH FRACTION ******************!
!!!!!!!!!!!!!!!!!! DEFINITION OF TOP AND BASE LAYER !!!!!!!!!!!!!!!!!!!!!!!!!
! low level < 3.2km
! 3.2 <= mid level < 6.5km
! high level >= 6.5km

!  680 hPa ===> 3.5km avec équilibre hydrostatique 
                        !  P=P0.exp(-z/H) , H=8.5
           !  440 hPa ===> 7.2km
          ! top lvl of isccp 




print *, 'diagnostic fraction nuage subgrid'


! Allocate / initialization of instantaneous isccp variables
 allocate(isccplow(it),isccpmid(it),isccphigh(it),colcloud(it))!, colclear(it))
 allocate(watercloud(altmax,it),icecloud(altmax,it),uncloud(altmax,it,catmax),phasecloud(altmax,it))
allocate(cftemp(tempmax-1,it),cftempliq(tempmax-1,it),cftempice(tempmax-1,it))


 allocate(height(it),height2(it),isccpliq(4,it),isccpice(4,it),isccpun(4,it,catmax))


   height(:)=0; height2(:)=0;
   icecloud(:,:)=0; watercloud(:,:)=0; uncloud(:,:,:)=0;phasecloud(:,:)=0;
   isccplow(:)=0; isccpmid(:)=0; isccphigh(:)=0;
   colcloud(:)=0; 
   isccpice(:,:)=0; isccpliq(:,:)=0; isccpun(:,:,:)=0;

 cftemp(:,:)=0; 
 cftempliq(:,:)=0; 
 cftempice(:,:)=0;


!! looking for limit before the 2 different SNR
do  iz=altmax,1,-1 
   if(altmod(iz).lt.8.16)then
      toplvlsat1=iz+1
      exit
   endif      
enddo

 



do i=1,it       !!!!! BEGIN OF IT LOOP 

   nanprof=0
   perptmp1 =0


! CLoud Phase diagnostic with the equation of the liquid relation between ATB
! and ATBr: ATBr = 1.3919 * ATB² + 0.0176 * ATB

 do  iz=altmax,toplvlsat1,-1  
   if(cloudfraction(iz,i).gt.0.)then
       if(perpmoy(iz,i).gt.parmoy(iz,i))then   ! noisy point = unphysical value

          uncloud(iz,i,5)=uncloud(iz,i,5)+1.   
          phasecloud(iz,i)=7.
          
       else

!! Curve for discrimination between ice & liquid particles
           perptmp1 = (pr2moy(iz,i)**5)*alpha50 + (pr2moy(iz,i)**4)*beta50 + (pr2moy(iz,i)**3)*gamma50 + (pr2moy(iz,i)**2)*delta50 + pr2moy(iz,i)*epsilon50 + zeta50

!! Curve for discrimination between ho particles & other
           perptmp2 = pr2moy(iz,i)*Ahoi + Bhoi

            if( (perpmoy(iz,i)-perptmp1).ge.0. )then
                 if (tempmoy(iz,i).gt.0 ) then
                    uncloud(iz,i,3)=uncloud(iz,i,3)+1.    !fake ice particles  
                    watercloud(iz,i)=watercloud(iz,i)+1.  !fake ice particles 
                    phasecloud(iz,i)=5.
                 else
                    icecloud(iz,i)=icecloud(iz,i)+1.      ! ice particles
                    phasecloud(iz,i)=2.
                 endif

            elseif( (perpmoy(iz,i)-perptmp2).ge.0. )then
                    if (tempmoy(iz,i).gt.-42 ) then
                       watercloud(iz,i)=watercloud(iz,i)+1.
                       phasecloud(iz,i)=1.
                    else
                       icecloud(iz,i)=icecloud(iz,i)+1.      ! fakeliq particles
                       uncloud(iz,i,2)=uncloud(iz,i,2)+1.    !fakeliq particles  
                       phasecloud(iz,i)=4.                   !fakeliq particles
                    endif
                 else
                   uncloud(iz,i,4)=uncloud(iz,i,4)+1. !! horizontally oriented particles
                   phasecloud(iz,i)=6.              !! horizontally oriented particles
            endif
        
         endif
    endif
  enddo 


 do  iz=toplvlsat1-1,1,-1  
   if(cloudfraction(iz,i).gt.0.)then
       if(perpmoy(iz,i).gt.parmoy(iz,i))then   ! noisy point = unphysical value

          uncloud(iz,i,5)=uncloud(iz,i,5)+1.  
          phasecloud(iz,i)=7.           
       else

!! Curve for discrimination between ice & liquid particles
           perptmp1 = (pr2moy(iz,i)**5)*alpha50 + (pr2moy(iz,i)**4)*beta50 + (pr2moy(iz,i)**3)*gamma50 + (pr2moy(iz,i)**2)*delta50 + pr2moy(iz,i)*epsilon50 + zeta50

!! Curve for discrimination between ho particles & other
           perptmp2 = pr2moy(iz,i)*Ahoi + Bhoi
   
            if( (perpmoy(iz,i)-perptmp1).ge.0. )then
                 if (tempmoy(iz,i).gt.0 ) then
                    uncloud(iz,i,3)=uncloud(iz,i,3)+1.    !fake ice particles
                    watercloud(iz,i)=watercloud(iz,i)+1.  !fake ice particles   
                    phasecloud(iz,i)=5.
                 else
                    icecloud(iz,i)=icecloud(iz,i)+1.
                    phasecloud(iz,i)=2.
                 endif

            elseif( (perpmoy(iz,i)-perptmp2).ge.0. )then
                    if (tempmoy(iz,i).gt.-42 ) then
                       watercloud(iz,i)=watercloud(iz,i)+1.
                       phasecloud(iz,i)=1.
                    else
                       icecloud(iz,i)=icecloud(iz,i)+1.      ! fakeliq particles
                       uncloud(iz,i,2)=uncloud(iz,i,2)+1.    !fake liq particles  
                       phasecloud(iz,i)=4.                   !fake liq particles
                    endif
                 else
                   uncloud(iz,i,4)=uncloud(iz,i,4)+1.  !! horizontally oriented particles
                   phasecloud(iz,i)=6.                 !! horizontally oriented particles
            endif

       toplvlsat2=0
    if( srmoy(iz,i).gt.30)then
       toplvlsat2=iz

       if(toplvlsat2.gt.altmax-2)then
          print *, toplvlsat2
       endif

       goto 99 
    endif

   endif
  endif
 enddo 

99 continue

if(toplvlsat2.ne.0)then

! Other level below cloud with SR>30 phase = undefined
do iz=toplvlsat2-1,1,-1
    if(cloudfraction(iz,i).gt.0.)then
      uncloud(iz,i,1)=uncloud(iz,i,1)+1.
      phasecloud(iz,i)=3.  !! undefined cloud
    endif
enddo
       toplvlsat2=0
endif



do  iz=1,altmax  
	if(cloudfraction(iz,i).gt.0.)then
		icewaterres=watercloud(iz,i)+icecloud(iz,i)+uncloud(iz,i,1)+uncloud(iz,i,4)+uncloud(iz,i,5)
		if((icewaterres.gt.1.).or.(icewaterres.eq.0.))then 
			print *, 'error sum phase=',icewaterres
			print *, watercloud(iz,i),icecloud(iz,i),uncloud(iz,i,:)
			!fixme
			stop
		endif
     
	endif
	enddo



altend=0
altstart=0

if(nol.eq.1)then


!!! NON OVERLAP MODE
!!! select only clouds in the highest isccp layer
!!! if isccpmid=1 then isccphigh must be 0
!!! A isccp layer can not have a cloud above 
B34: do  iz=altmax,1,-1  
     if(cloudfraction(iz,i).gt.0.)then


! Search the isscp high cloud fraction
         if (iz.ge.topmidl) then
     altstart=altmax
     altend=topmidl
            exit B34
         endif

! Search the isscp mid cloud fraction
         if ((iz.ge.toplowl).and.(iz.lt.topmidl)) then
         altstart=topmidl-1
         altend=toplowl
            exit B34
         endif

! Search the isscp low cloud fraction
        if (iz.lt.toplowl) then
            altstart=toplowl-1
            altend=1
            exit B34
       endif

    endif

enddo B34

else



altend=1
altstart=altmax

endif

if(altend.ne.0)then   ! to avoid clear sky profil error
indbase=0
 do  iz=altend,altstart  

   if(cloudfraction(iz,i).gt.0.)then

! Calculation of cloudfraction as a function of the temperature
  do itemp=1,tempmax-1
    if ( (tempmoy(iz,i).gt.tempmod(itemp)).and.     &
       (tempmoy(iz,i) .le.tempmod(itemp+1)) )then

       cftemp(itemp,i)=cftemp(itemp,i)+1   ! cf tot

       if( watercloud(iz,i).ne.0 )then   
          cftempliq(itemp,i)=cftempliq(itemp,i)+1  ! cf liquid
       elseif( icecloud(iz,i).ne.0 )then  
          cftempice(itemp,i)=cftempice(itemp,i)+1  ! cf ice
       endif

    endif
  enddo


! Indice de depolarisation (indices: 0=clear-sky, 1=liquid only; 2=mixte; 3=ice only)
       
       depoltmp=0
       icewaterres=0

! calculation of the depol of the straight line  Depol = 0.0028*SR + 0.0458
! rajouter un critere de temperature pour suprimer les fake liquide ?
! rajouter critere d'altitude ?

! CLoud Phase diagnostic with the equation of the liquid relation between ATB
! and ATBr: ATBr = 1.3919 * ATB² + 0.0176 * ATB
!  




! Search the isscp low cloud fraction
   if (iz.lt.toplowl) then
      isccplow(i)=1

      if( icecloud(iz,i).gt.0 )then
         isccpice(1,i)=1 
      endif

      if( watercloud(iz,i).gt.0 )then
         isccpliq(1,i)=1
      endif
      
      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(1,i,icat)=1
         endif
      enddo

   endif

! Search the isscp mid cloud fraction
   if ((iz.ge.toplowl).and.(iz.lt.topmidl)) then
      isccpmid(i)=1
            
      if(icecloud(iz,i).gt.0)then
         isccpice(2,i)=1 
      endif

      if(watercloud(iz,i).gt.0)then
         isccpliq(2,i)=1
      endif

      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(2,i,icat)=1
         endif
      enddo

   endif

! Search the isscp high cloud fraction
   if (iz.ge.topmidl) then
      isccphigh(i)=1

height(i)=altmid(iz)+0.24

indbase=indbase+1
if(indbase.eq.1)then
  height2(i)=altmid(iz)-0.24
endif

      if(icecloud(iz,i).gt.0)then
         isccpice(3,i)=1 
      endif
      if(watercloud(iz,i).gt.0)then
         isccpliq(3,i)=1
      endif

      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(3,i,icat)=1
         endif
      enddo

   endif

! Search cloud fraction on the column
            colcloud(i)=1

      if(icecloud(iz,i).gt.0)then
         isccpice(4,i)=1 
      endif

      if(watercloud(iz,i).gt.0)then
         isccpliq(4,i)=1
      endif

      do icat=1,catmax
         if(uncloud(iz,i,icat).gt.0)then
            isccpun(4,i,icat)=1
         endif
      enddo

    endif  !! endif cloudfraction > 0

         if((nanfraction(iz,i).eq.1.).or.(rejfraction(iz,i).eq.1.).or.(sefraction(iz,i).eq.1.))then
            nanprof=nanprof+1
         endif



  enddo
endif 


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                  DAILY AVERAGE OF CLOUD DIAGNOSTICS                      *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


       
!******************* READING THE LAT LON MODEL GRID FILES *******************!

!---------------------------- lat lon grid file -----------------------------!
!                                                                            !
!   There is 3 grid in this program :                                        !
!                                                                            !
!       -CFMIP : 1° x 1° x 41levels each 480m from 0 to 19.2km               !
!                (361,181,41)                                                !
!       -LMDZ  : 3.75° x 2.53° x 19levels from 0 to 40.4823km                !
!                (96,72,19)                                                  !
!       -NASA  : 5° x 5° x 41levels each 480m from 0 to 19.2km               !
!               (73,37,41)                                                   !
!                                                                            !
!----------------------------------------------------------------------------!


!************************** BEGIN THE DAILY AVERAGE *************************!

!---------------------------- lat lon grid file -----------------------------!
!                                                                            !
!   increase the daily variables over the lat/lon (& alt & diag) grid        !
!                                                                            !
!----------------------------------------------------------------------------!

! First level always equal to 0 because the average is perfomed between the 
! level i and i-1.



  do ilon=1,lonmax-1  !longitude
     if( (lon(i).ge.lonmod(ilon)) .and. (lon(i).lt.lonmod(ilon+1)) )then
         
    do ilat=1,latmax-1  !latitude
        if ( (lat(i).ge.latmod(ilat)) .and. (lat(i).lt.latmod(ilat+1)) )then

            indtot(ilon,ilat)=indtot(ilon,ilat)+1
             
             nanprof=0
             do ialt=1,altmax
                if((nanfraction(ialt,i).eq.1.).or.(sefraction(ialt,i).eq.1.) &
                   .or.(rejfraction(ialt,i).eq.1.))then
                   nanprof=nanprof+1
                endif
             enddo
             
             nansat=0  
             nanlow=0
             do ialt=1,toplowl-1
                if(satfraction(ialt,i).eq.1.)then
                   nansat=nansat+1       ! fully att layer
                endif
                if(sefraction(ialt,i).eq.1.)then
                   nanlow=nanlow+1       ! under SE layer
                endif
             enddo
         
             nanmid=0
             do ialt=toplowl,topmidl-1
                if(satfraction(ialt,i).eq.1.)then
                   nanmid=nanmid+1
                endif
             enddo

          ! Cloudtop Height average (GEWEX option)
          heightday(ilat,ilon,jour)=heightday(ilat,ilon,jour)+height(i)
          heightday2(ilat,ilon,jour)=heightday2(ilat,ilon,jour)+height2(i)
          if( height(i).gt.0)then
             indheight(ilat,ilon,jour)=indheight(ilat,ilon,jour)+1
          endif

do nphase = 1,4
isccpliqday(ilat,ilon,jour,nphase)=isccpliqday(ilat,ilon,jour,nphase)+isccpliq(nphase,i)
isccpiceday(ilat,ilon,jour,nphase)=isccpiceday(ilat,ilon,jour,nphase)+isccpice(nphase,i)
  do icat=1,catmax
isccpunday(ilat,ilon,jour,nphase,icat)=isccpunday(ilat,ilon,jour,nphase,icat)+isccpun(nphase,i,icat)
  enddo
enddo

!!!! A RAJOUTER LES ISCCPblablaDAY
!!!! ENREGISTRER LES 2 TYPE DE ISCCP en occurrence ?

         ! colclearday(ilat,ilon,jour)=colclearday(ilat,ilon,jour)+colclear(i)
          
          ! select only the profil .ne. NaN & Low layer not fully attenuated
          if( (nanprof.le.(altmax-4)).and.(nansat.lt.(toplowl-1)).and. &
             (nanlow.lt.(toplowl-1)) )then
          isccplowday(ilat,ilon,jour)=isccplowday(ilat,ilon,jour)+isccplow(i)
          isccpinddaylow(ilat,ilon,jour)=isccpinddaylow(ilat,ilon,jour)+1
          endif

          ! select only the profil .ne. NaN & Mid layer not fully attenuated
          if((nanprof.le.(altmax-4)).and.(nanmid.lt.(topmidl-toplowl)))then
          isccpinddaymid(ilat,ilon,jour)=isccpinddaymid(ilat,ilon,jour)+1
          isccpmidday(ilat,ilon,jour)=isccpmidday(ilat,ilon,jour)+isccpmid(i)
          endif

          ! select only the profil .ne. NaN
          if(nanprof.le.(altmax-4))then
          isccpindday(ilat,ilon,jour)=isccpindday(ilat,ilon,jour)+1  
          isccphighday(ilat,ilon,jour)=isccphighday(ilat,ilon,jour)+         &
                                       isccphigh(i)
          colcloudday(ilat,ilon,jour)=colcloudday(ilat,ilon,jour)+colcloud(i)    
          endif

          if( (isccplowday(ilat,ilon,jour)/isccpinddaylow(ilat,ilon,jour)) .gt. 1.00001)then
             print *, file6
             print *, i,lon(i),lat(i),nansat
             print *, isccplowday(ilat,ilon,jour),isccpinddaylow(ilat,ilon,jour)
          endif

 
 do  iz=1,altmax  
     if(cloudfraction(iz,i).gt.0.)then

if(indicetemp(iz,i).gt.0)then
       coltemp(ilat,ilon,jour)=coltemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
       indcoltemp(ilat,ilon,jour)=indcoltemp(ilat,ilon,jour)+1

! Search the isscp low cloud fraction
        if (iz.lt.toplowl) then
            lowtemp(ilat,ilon,jour)=lowtemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
            indlowtemp(ilat,ilon,jour)=indlowtemp(ilat,ilon,jour)+1
       endif

! Search the isscp mid cloud fraction
         if ((iz.ge.toplowl).and.(iz.lt.topmidl)) then
            midtemp(ilat,ilon,jour)=midtemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
            indmidtemp(ilat,ilon,jour)=indmidtemp(ilat,ilon,jour)+1
         endif

! Search the isscp high cloud fraction
         if (iz.ge.topmidl) then
             hightemp(ilat,ilon,jour)=hightemp(ilat,ilon,jour)+(tempmoy(iz,i)/indicetemp(iz,i))
            indhightemp(ilat,ilon,jour)=indhightemp(ilat,ilon,jour)+1
         endif
endif

      endif
 enddo


! cloudfraction day as a function of temperature
  do itemp=1,tempmax-1
     cftempday(ilat,ilon,itemp,jour)=cftempday(ilat,ilon,itemp,jour)+cftemp(itemp,i)
     cftempliqday(ilat,ilon,itemp,jour)=cftempliqday(ilat,ilon,itemp,jour)+   &
                                   cftempliq(itemp,i)
     cftempiceday(ilat,ilon,itemp,jour)=cftempiceday(ilat,ilon,itemp,jour)+   &
                                   cftempice(itemp,i)
  enddo


            do ialt=altmax,1,-1  
            indnan(ilat,ilon,ialt)=indnan(ilat,ilon,ialt)+1
            cloudfractday(ilat,ilon,ialt,jour)=                              &
            cloudfractday(ilat,ilon,ialt,jour)+cloudfraction(ialt,i)
            clearfractday(ilat,ilon,ialt,jour)=                              &
            clearfractday(ilat,ilon,ialt,jour)+clearfraction(ialt,i)
            uncertfractday(ilat,ilon,ialt,jour)=                             &
            uncertfractday(ilat,ilon,ialt,jour)+uncertfraction(ialt,i)
            icecloudfractday(ilat,ilon,ialt,jour)=                            &
            icecloudfractday(ilat,ilon,ialt,jour)+icecloud(ialt,i)
            watercloudfractday(ilat,ilon,ialt,jour)=                          &
            watercloudfractday(ilat,ilon,ialt,jour)+watercloud(ialt,i)

            do icat=1,catmax
            uncloudfractday(ilat,ilon,ialt,jour,icat)=                       &
            uncloudfractday(ilat,ilon,ialt,jour,icat)+uncloud(ialt,i,icat)
            enddo

            if((nanfraction(ialt,i).ne.1).and.(sefraction(ialt,i).ne.1).and.(satfraction(ialt,i).ne.1).and.(rejfraction(ialt,i).ne.1))then
            indday(ilat,ilon,ialt,jour)=indday(ilat,ilon,ialt,jour)+1

  do itemp=1,tempmax-1
               if ( (tempmoy(ialt,i).gt.tempmod(itemp)).and.     &
                   (tempmoy(ialt,i) .le.tempmod(itemp+1)) )then
                   indcftemp(ilat,ilon,itemp,jour)=indcftemp(ilat,ilon,itemp,jour)+1
                endif
  enddo
            ! indice of cloudfraction day as a function of temperature

            if(indday(ilat,ilon,ialt,jour).lt.sum(uncloudfractday(ilat,ilon,ialt,jour,:)))then
               print *, 'error indice < un_cloud',i,ialt
               ! FIXME
               stop
            endif
        
            endif

            if(srmoy(ialt,i).eq.srmod(1))then
                diagSR(ilon,ilat,ialt,1)=diagSR(ilon,ilat,ialt,1)+1                
            endif
            if(srmoy(ialt,i).eq.srmod(2))then
                diagSR(ilon,ilat,ialt,2)=diagSR(ilon,ilat,ialt,2)+1
            endif


                  do idiag=3,diagmax-1
                     if ( (srmoy(ialt,i).ge.srmod(idiag)).and.             &
                        (srmoy(ialt,i).lt.srmod(idiag+1)) )then
                        diagSR(ilon,ilat,ialt,idiag)=                   &
                        diagSR(ilon,ilat,ialt,idiag)+1
                     endif
                  enddo

                  do idiag=8,diagmax-1
                     if ( (srmoy(ialt,i).ge.srmod(idiag)).and.             &
                        (srmoy(ialt,i).lt.srmod(idiag+1)) )then
                       if(icecloud(ialt,i).ne.0)then
                        diagSRpha(ilon,ilat,ialt,idiag-7,2)=                   &
                        diagSRpha(ilon,ilat,ialt,idiag-7,2)+1

                        elseif(watercloud(ialt,i).ne.0)then
                        diagSRpha(ilon,ilat,ialt,idiag-7,1)=                   &
                        diagSRpha(ilon,ilat,ialt,idiag-7,1)+1

                        elseif(uncloud(ialt,i,1).ne.0)then
                        diagSRpha(ilon,ilat,ialt,idiag-7,3)=                   &
                        diagSRpha(ilon,ilat,ialt,idiag-7,3)+1
   
                        endif
                     endif
                  enddo


              enddo
            endif
         enddo
      endif
   enddo



enddo !!!!!!!!!! END IT LOOP



!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                   RECORD INSTANTANEOUS PHASE FILE                        *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!
if('instant_switch'.eq.'on')then
 call instant_phase(file4,altmax,it,phasecloud)
endif


do ilon=1,lonmax-1  !latitude
        do ilat=1,latmax-1  !longitude
  
          if(indtot(ilon,ilat).gt.0)then
           indtotmean(ilon,ilat)=indtotmean(ilon,ilat)+1
          endif
        enddo
enddo

indtot(:,:)=0;


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                           WRITE OUTPUT FILES                             *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!----------------------------- Save output file -----------------------------!
!                                                                            !
! Sauve = "wrf" when model WRF is used, then the data are saved in ASCII       !
! format, and for every hdf calipso file. Only the matching WRF grid values  !
! are saved.                                                                 !
!     - Lat : 199values between 15° to 55°                                   !
!     - Lon : 179values between -18° to 36°                                  !
!     - alt : 48levels from 0 to 35.5km                                      !
!                                                                            !
! Sauve = "chim" when model CHIMERE is used, then the data are saved in ASCII  !
! format, and for every hdf calipso file. Only the matching CHIMERE grid     !
! values are saved.                                                          !
!     - Lat : 71values between -5° to 65° (every 1°)                         !
!     - Lon : 181values between -100° to 80° (every 1°)                      !
!     - alt : 281levels from 0 to 14km (every 50m)                           !
!                                                                            !          
! Sauve = "lmdz" when LMDZ is used, the data are saved in a single netcdf file !
! when the last hdf calipso file of the period have been read                !        
!                                                                            !
!----------------------------------------------------------------------------!
622 continue

sauve=model

SELECT CASE (sauve)


!**************************** WRF OUTPUT FORMAT *****************************! 
 CASE ("wrf")

comptpf=0
do i=1,it
     if ((lat(i).gt.15.0).and.(lat(i).lt.55.0))then
           if((lon(i).gt.-18.0).and.(lon(i).lt.36.0))then
              comptpf=comptpf+1
           endif
      endif
enddo

allocate(latwrf(comptpf),lonwrf(comptpf),SRwrf(altmax,comptpf),DEPOLwrf(altmax,comptpf),CRwrf(altmax,comptpf),SEwrf(comptpf),temps2wrf(comptpf))

if(comptpf.eq.0)then
   goto 624  !track out from the domain
else

comptpf=0
do i=1,it
     if ((lat(i).gt.15.0).and.(lat(i).lt.55.0))then
           if((lon(i).gt.-18.0).and.(lon(i).lt.36.0))then
              comptpf=comptpf+1
              SEwrf(comptpf)=SE(i)
              temps2wrf(comptpf)=temps2(i)
              latwrf(comptpf)=lat(i)
              lonwrf(comptpf)=lon(i)
              do j=1,altmax
                 SRwrf(j,comptpf)=srmoy(j,i)
                 CRwrf(j,comptpf)=crmoy(j,i)
                 DEPOLwrf(j,comptpf)=depolmoy(j,i)
              enddo
           endif
      endif
enddo


  write(numfichc,'(i4)')numfich
  write(datec,'(I5)')date

! Name of instantaneous file, one by hdf calipso file  
  file4='SR_CR_DEPOL_200'//datec//'_'//trim(switch)//'_'//trim(gcm)//'_'//    &
  trim(ADJUSTL(numfichc))//'.nc'

  call SR_CR_DR_2nc(file4,altmid,altmod,resd,altmax,switch,gcm,comptpf,latwrf,lonwrf,SEwrf,temps2wrf,&
                  SRwrf,CRwrf,DEPOLwrf)

endif

624 continue



! Deallocate SDS/META variables
deallocate(atb,stat = OK_buffer)
deallocate(atb2,stat = OK_buffer)
deallocate(perp,stat = OK_buffer)
deallocate(mol,stat = OK_buffer)
deallocate(pres,stat = OK_buffer) 
deallocate(SE,SEwrf,stat = OK_buffer)
deallocate(mol2,stat = OK_buffer)
deallocate(mol3,stat = OK_buffer)
deallocate(mol4)
deallocate(temp,stat = OK_buffer)
deallocate(temp2,stat = OK_buffer)

if(alt_pres=='pressure')then
   deallocate(pres2,stat = OK_buffer)
endif
deallocate(altl,stat = OK_buffer)
deallocate(altm,stat = OK_buffer)
deallocate(temps2,temps2wrf,stat = OK_buffer)
deallocate(lonwrf,latwrf,SRwrf,CRwrf,DEPOLwrf)


! end of select "wrf"

!continue 


!************************** CHIMERE OUTPUT FORMAT ***************************!
CASE ("chim")

print *, 'Regrid of data Done'
print *, ''
print *, 'Save the data'

do k=1,lonmax
   do j=1,latmax 
      if(mheure(j,k).ne.0)then
      box=box+1
      if(lonmod(k)==12)then
       
endif

      write(11,'(4(2x,I6),2(2x,F10.2),(2x,E13.6))'),numfich,box,month,day,   &
      lonmod(k),latmod(j),mheure(j,k)/indiceh(j,k)
         do iz=1,altmax
      
            enddo
      endif
      
   enddo
enddo


!**************************** LMDZ OUTPUT FORMAT ****************************!
CASE ("lmdz")

  print *,'model==lmdz, daily files creation skipped.'
  goto 666

CASE DEFAULT

  ! model is not in [lmdz, chimere, wrf] - should never happen

  print *,'invalid select : model = ', sauve
  print *,'whatever, chugging along...'

END SELECT

666 continue
647 continue


!***** DEALLOCATE SDS OUTPUT AND INSTANT LMDZ OUTPUT VAR AND CLOSE FILES ****!

deallocate(pr2moy,molmoy, stat = OK_buffer)!
deallocate(srmoy)
deallocate(lat,stat = OK_buffer)
deallocate(lon,stat = OK_buffer)
deallocate(indice,indicem,indice2)
deallocate(indicep,indicep2,depolmoy,parmoy,perpmoy,pr2moy2,crmoy)
deallocate(tempmoy,indicetemp)

if(model=='wrf')then
 deallocate(indice,indicem,indice2)
deallocate(indicep,indicep2,crmoy,depolmoy,parmoy,perpmoy,pr2moy2)

endif


if(model=='chimere')then
   deallocate(latmod,lonmod,altmod,stat = OK_buffer)
   deallocate(mheure,stat = OK_buffer)
   deallocate(indiceh,stat = OK_buffer)
   deallocate(indice,indicem,indice2,stat = OK_buffer)
endif

if(model=='lmdz')then

if(alt_pres=='pressure')then
 deallocate(SEp, stat = OK_buffer)
endif
 deallocate(isccplow,isccpmid,isccphigh,colcloud)
 deallocate(cloudfraction,clearfraction,satfraction,uncertfraction,          &
            nanfraction,sefraction,rejfraction)
 deallocate(icecloud,watercloud,uncloud,phasecloud)
 deallocate(isccpliq,isccpice,isccpun)
 deallocate(cftemp,cftempliq,cftempice)
 deallocate(height,height2)
endif

print *, 'Deallocate buffers done'


887 continue
print *, ''

print *, 'go to the next file'


  ! Temps elapsed final
  call system_clock(count=t2, count_rate=ir)
  tempstot=real(t2 - t1,kind=4)/real(ir,kind=4)

  ! Temps CPU de calcul final
  call cpu_time(t_cpu_1)
  t_cpu = t_cpu_1 - t_cpu_0

         
 print '(//,3X,"Temps elapsed       : ",1PE10.3," sec.",/, &
           & 3X,"Temps CPU           : ",1PE10.3," sec.",/,//)', &
           tempstot,t_cpu




goto 888          ! next Calipso file


999 continue 
 
! empty file list checking
if(file2(1:1).gt.' ')then
continue
else
print *, 'hdf file list empty'
print *, file2(1:1)
stop
endif


!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                  MONTHLY AVERAGE CLOUD DIAGNOSTICS                       *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!----------------------------- monthly average ------------------------------!
!                                                                            !
!  If there is no more hdf CALIPSO file to read, then the monthly average    !
!  can begin. This module is operating only in the lmdz case, that is to say !
!  when we have to save the monthly diagnostic.                              !
!  In the other case (chimere & wrf), the program ends here.                 !
!                                                                            !
!----------------------------------------------------------------------------!

SELECT CASE (model)

 CASE ("lmdz")

!****************************************************************************!
!*!!!!!!!!!!!!!! PART I : CLOUDY LOW MID HIGH MAP FILES !!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


! Allocation / initialization of isscp monthly variables
print *, 'allocation / initialization of isccp monthly variables'

   allocate(monthheight(lonmax-1,latmax-1),indmonthheight(lonmax-1,latmax-1),      &
            monthheight2(lonmax-1,latmax-1))
   allocate(monthisccplow(lonmax-1,latmax-1),monthisccpmid(lonmax-1,latmax-1),       &
            monthisccphigh(lonmax-1,latmax-1))
   allocate(monthcolcloud(lonmax-1,latmax-1), monthcolclear(lonmax-1,latmax-1),       &
            isccpdaypermonthlow(lonmax-1,latmax-1),                                  &
            isccpdaypermonthmid(lonmax-1,latmax-1),                                  &           
            isccpdaypermonth(lonmax-1,latmax-1))
   allocate(hlow(lonmax-1,latmax-1,histmax-1),hmid(lonmax-1,latmax-1,histmax-1), &
            hhigh(lonmax-1,latmax-1,histmax-1),hcol(lonmax-1,latmax-1,histmax-1),&
            hheight(lonmax-1,latmax-1,histmax2-1))

   allocate(hlowtemp(lonmax-1,latmax-1,histtempmax2),hmidtemp(lonmax-1,latmax-1,histtempmax2), &
            hhightemp(lonmax-1,latmax-1,histtempmax2),hcoltemp(lonmax-1,latmax-1,histtempmax2))
   allocate(monthlowtemp(lonmax-1,latmax-1),indmonthlowtemp(lonmax-1,latmax-1), &
            monthmidtemp(lonmax-1,latmax-1),indmonthmidtemp(lonmax-1,latmax-1), &
            monthhightemp(lonmax-1,latmax-1),indmonthhightemp(lonmax-1,latmax-1), &
            monthcoltemp(lonmax-1,latmax-1),indmonthcoltemp(lonmax-1,latmax-1))


allocate(monthisccpliq(lonmax-1,latmax-1,4),monthisccpice(lonmax-1,latmax-1,4))
allocate(monthisccpun(lonmax-1,latmax-1,4,catmax))
allocate(monthisccpphase(lonmax-1,latmax-1,4))
allocate(indmonthphase(lonmax-1,latmax-1,4),indmonthphase2(lonmax-1,latmax-1,4))
allocate(inddayphase(latmax-1,lonmax-1,daymax,4))

inddayphase(:,:,:,:)=0;
monthisccpliq(:,:,:)=0; monthisccpice(:,:,:)=0; 
monthisccpun(:,:,:,:)=0; 
monthisccpphase(:,:,:)=0;indmonthphase(:,:,:)=0;indmonthphase2(:,:,:)=0

histmod(:)=0; histmod2(:)=0;
hlow(:,:,:)=0;hmid(:,:,:)=0;hhigh(:,:,:)=0;hcol(:,:,:)=0;hheight(:,:,:)=0;
monthheight(:,:)=0; monthheight2(:,:)=0; indmonthheight(:,:)=0;
monthisccplow(:,:)=0;monthisccpmid(:,:)=0;monthisccphigh(:,:)=0;
monthcolcloud(:,:)=0;isccpdaypermonth(:,:)=0;monthcolclear(:,:)=0
isccpdaypermonthlow(:,:)=0;isccpdaypermonthmid(:,:)=0;
histtempmod(:)=0; histtempmod2(:)=0;

hlowtemp(:,:,:)=0;monthlowtemp(:,:)=0;indmonthlowtemp(:,:)=0;
hmidtemp(:,:,:)=0;monthmidtemp(:,:)=0;indmonthmidtemp(:,:)=0;
hhightemp(:,:,:)=0;monthhightemp(:,:)=0;indmonthhightemp(:,:)=0;
hcoltemp(:,:,:)=0;monthcoltemp(:,:)=0;indmonthcoltemp(:,:)=0;

histtempmod(1)=150
histtempmod(2)=180
histtempmod(histtempmax)=320

do ihisttemp=2,histtempmax-2
   histtempmod(ihisttemp+1)=histtempmod(ihisttemp)+5
enddo


do ihisttemp=1,histtempmax2
   histtempmod2(ihisttemp)=(histtempmod(ihisttemp)+histtempmod(ihisttemp+1))/2
enddo


do ihist=1,histmax-2
   histmod(ihist+1)=histmod(ihist)+0.1
enddo
 
histmod(histmax)=1.01

do ihist=1,histmax2-1
   histmod2(ihist+1)=histmod2(ihist)+1
enddo




!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!
! changement à rajouter !
! normaliser le calcul de la phase par jour


!********* CALCULATION OF DAILY DIAGNOSTIC WITH MATCHING INDEXES ************!
! moyenne journalière effectuée avant la moyenne mensuelle pour donner le même
! poid à tous les jours malgré la quantité variable de donnée d'un jour à 
! l'autre.

      do jour=1,31 
         do ilon=1,lonmax-1
           do ilat=1,latmax-1  

             if ( indlowtemp(ilat,ilon,jour).gt.0 ) then
           lowtemp(ilat,ilon,jour)=                                         &
           (lowtemp(ilat,ilon,jour)/indlowtemp(ilat,ilon,jour))+273.15
             endif
            if ( indmidtemp(ilat,ilon,jour).gt.0 ) then
           midtemp(ilat,ilon,jour)=                                         &
           (midtemp(ilat,ilon,jour)/indmidtemp(ilat,ilon,jour))+273.15
             endif
            if ( indhightemp(ilat,ilon,jour).gt.0 ) then
           hightemp(ilat,ilon,jour)=                                         &
           (hightemp(ilat,ilon,jour)/indhightemp(ilat,ilon,jour))+273.15
             endif
            if ( indcoltemp(ilat,ilon,jour).gt.0 ) then
           coltemp(ilat,ilon,jour)=                                         &
           (coltemp(ilat,ilon,jour)/indcoltemp(ilat,ilon,jour))+273.15
             endif


             if ( isccpindday(ilat,ilon,jour).gt.0 ) then
           isccphighday(ilat,ilon,jour)=                                     &
           isccphighday(ilat,ilon,jour)/isccpindday(ilat,ilon,jour)
           colcloudday(ilat,ilon,jour)=                                      &
           colcloudday(ilat,ilon,jour)/isccpindday(ilat,ilon,jour)

if (isccpiceday(ilat,ilon,jour,3)+isccpliqday(ilat,ilon,jour,3)+         &
    +sum(isccpunday(ilat,ilon,jour,3,4:5))+isccpunday(ilat,ilon,jour,3,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,3)+isccpliqday(ilat,ilon,jour,3)+          &
          isccpunday(ilat,ilon,jour,3,1)+ sum(isccpunday(ilat,ilon,jour,3,4:5))

          isccpliqday(ilat,ilon,jour,3)= isccpliqday(ilat,ilon,jour,3)/            &
          isccptemp * isccphighday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,3)= isccpiceday(ilat,ilon,jour,3)/            &
          isccptemp * isccphighday(ilat,ilon,jour)
do icat=1,catmax
          isccpunday(ilat,ilon,jour,3,icat)= isccpunday(ilat,ilon,jour,3,icat)/            &
          isccptemp * isccphighday(ilat,ilon,jour)
enddo

endif 

if (isccpiceday(ilat,ilon,jour,4)+isccpliqday(ilat,ilon,jour,4)+         &
    +sum(isccpunday(ilat,ilon,jour,4,4:5))+isccpunday(ilat,ilon,jour,4,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,4)+isccpliqday(ilat,ilon,jour,4)+          &
          isccpunday(ilat,ilon,jour,4,1)+ sum(isccpunday(ilat,ilon,jour,4,4:5))

          isccpliqday(ilat,ilon,jour,4)= isccpliqday(ilat,ilon,jour,4)/            &
          isccptemp * colcloudday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,4)= isccpiceday(ilat,ilon,jour,4)/            &
          isccptemp * colcloudday(ilat,ilon,jour)

do icat=1,catmax
          isccpunday(ilat,ilon,jour,4,icat)= isccpunday(ilat,ilon,jour,4,icat)/      &
          isccptemp * isccphighday(ilat,ilon,jour)
enddo

endif

             endif

             if ( indheight(ilat,ilon,jour).gt.0 ) then
                heightday(ilat,ilon,jour)=                                   &
                heightday(ilat,ilon,jour)/indheight(ilat,ilon,jour)                      
                heightday2(ilat,ilon,jour)=                                   &
                heightday2(ilat,ilon,jour)/indheight(ilat,ilon,jour)                  
             endif

             if ( isccpinddaylow(ilat,ilon,jour).gt.0 ) then
           isccplowday(ilat,ilon,jour)=                                      &
           isccplowday(ilat,ilon,jour)/isccpinddaylow(ilat,ilon,jour)

if (isccpiceday(ilat,ilon,jour,1)+isccpliqday(ilat,ilon,jour,1)+         &
    +sum(isccpunday(ilat,ilon,jour,1,4:5))+isccpunday(ilat,ilon,jour,1,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,1)+isccpliqday(ilat,ilon,jour,1)+          &
          isccpunday(ilat,ilon,jour,1,1)+ sum(isccpunday(ilat,ilon,jour,1,4:5))

          isccpliqday(ilat,ilon,jour,1)= isccpliqday(ilat,ilon,jour,1)/            &
          isccptemp * isccplowday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,1)= isccpiceday(ilat,ilon,jour,1)/            &
          isccptemp * isccplowday(ilat,ilon,jour)

do icat=1,catmax
          isccpunday(ilat,ilon,jour,1,icat)= isccpunday(ilat,ilon,jour,1,icat)/      &
          isccptemp * isccplowday(ilat,ilon,jour)
enddo

endif


             endif
             if ( isccpinddaymid(ilat,ilon,jour).gt.0 ) then
           isccpmidday(ilat,ilon,jour)=                                      &
           isccpmidday(ilat,ilon,jour)/isccpinddaymid(ilat,ilon,jour)

if (isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)+         &
    +sum(isccpunday(ilat,ilon,jour,2,4:5))+isccpunday(ilat,ilon,jour,2,1).gt.0)then
isccptemp=isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)+          &
          isccpunday(ilat,ilon,jour,2,1)+ sum(isccpunday(ilat,ilon,jour,2,4:5))

          isccpliqday(ilat,ilon,jour,2)= isccpliqday(ilat,ilon,jour,2)/            &
          isccptemp * isccpmidday(ilat,ilon,jour)
          isccpiceday(ilat,ilon,jour,2)= isccpiceday(ilat,ilon,jour,2)/            &
          isccptemp * isccpmidday(ilat,ilon,jour)

do icat=1,catmax
          isccpunday(ilat,ilon,jour,2,icat)= isccpunday(ilat,ilon,jour,2,icat)/     &
          isccptemp * isccpmidday(ilat,ilon,jour)
enddo 

endif



if ((isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)).gt.isccpmidday(ilat,ilon,jour))then
print *, isccpmidday(ilat,ilon,jour),isccpiceday(ilat,ilon,jour,2)+isccpliqday(ilat,ilon,jour,2)
endif

             endif
             if ( isccplowday(ilat,ilon,jour) .gt. colcloudday(ilat,ilon,jour) )then
                print *, ilon,ilat,jour
                print *, isccplowday(ilat,ilon,jour),colcloudday(ilat,ilon,jour)
             endif


do nphase=1,4
   isccptemp=isccpiceday(ilat,ilon,jour,nphase)+isccpliqday(ilat,ilon,jour,nphase)
   if(isccptemp.gt.0.)then
      inddayphase(ilat,ilon,jour,nphase)=inddayphase(ilat,ilon,jour,nphase)+1.
      isccpphaseday(ilat,ilon,jour,nphase)=isccpiceday(ilat,ilon,jour,nphase)/isccptemp
   endif
enddo

           enddo
         enddo
      enddo


!************* INCREMENTING THE MONTHLY DIAGNOSTIC AND INDEXES **************!

      do jour=1,31 
         do ilon=1,lonmax-1
           do ilat=1,latmax-1  

            if ( indlowtemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (lowtemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (lowtemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hlowtemp(ilon,ilat,ihisttemp)=                                         &
              hlowtemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthlowtemp(ilon,ilat)=monthlowtemp(ilon,ilat)+                  &
                                     lowtemp(ilat,ilon,jour)
            indmonthlowtemp(ilon,ilat)=indmonthlowtemp(ilon,ilat)+1            
            endif

            if ( indmidtemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (midtemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (midtemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hmidtemp(ilon,ilat,ihisttemp)=                                         &
              hmidtemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthmidtemp(ilon,ilat)=monthmidtemp(ilon,ilat)+                  &
                                     midtemp(ilat,ilon,jour)
            indmonthmidtemp(ilon,ilat)=indmonthmidtemp(ilon,ilat)+1            
            endif

            if ( indhightemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (hightemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (hightemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hhightemp(ilon,ilat,ihisttemp)=                                         &
              hhightemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthhightemp(ilon,ilat)=monthhightemp(ilon,ilat)+                  &
                                     hightemp(ilat,ilon,jour)
            indmonthhightemp(ilon,ilat)=indmonthhightemp(ilon,ilat)+1            
            endif

            if ( indcoltemp(ilat,ilon,jour).gt.0 ) then
            do ihisttemp=1,histtempmax-1
              if( (coltemp(ilat,ilon,jour).ge.histtempmod(ihisttemp)) .and.  &
                (coltemp(ilat,ilon,jour).lt.histtempmod(ihisttemp+1)) ) then
              hcoltemp(ilon,ilat,ihisttemp)=                                         &
              hcoltemp(ilon,ilat,ihisttemp)+1
              endif
            enddo
            monthcoltemp(ilon,ilat)=monthcoltemp(ilon,ilat)+                  &
                                     coltemp(ilat,ilon,jour)
            indmonthcoltemp(ilon,ilat)=indmonthcoltemp(ilon,ilat)+1            
            endif


          if ( indheight(ilat,ilon,jour).gt.0 ) then
             do ihist=1,histmax2-1
                if( (heightday(ilat,ilon,jour).ge.histmod2(ihist)) .and.&
                    (heightday(ilat,ilon,jour).lt.histmod2(ihist+1)) ) then
                    hheight(ilon,ilat,ihist)=hheight(ilon,ilat,ihist)+1
                 endif
              enddo
             monthheight(ilon,ilat)=monthheight(ilon,ilat)+                  &
                                    heightday(ilat,ilon,jour)
             monthheight2(ilon,ilat)=monthheight2(ilon,ilat)+                  &
                                    heightday2(ilat,ilon,jour)
             indmonthheight(ilon,ilat)=indmonthheight(ilon,ilat)+1
          endif

          if ( isccpindday(ilat,ilon,jour).gt.0 ) then

              do ihist=1,histmax-1
                 if( (isccphighday(ilat,ilon,jour).ge.histmod(ihist)) .and.  &
                     (isccphighday(ilat,ilon,jour).lt.histmod(ihist+1)) ) then
                     hhigh(ilon,ilat,ihist)=hhigh(ilon,ilat,ihist)+1
                 endif

                 if( (colcloudday(ilat,ilon,jour).ge.histmod(ihist)) .and.  &
                     (colcloudday(ilat,ilon,jour).lt.histmod(ihist+1)) )then
                     hcol(ilon,ilat,ihist)=hcol(ilon,ilat,ihist)+1
                  endif
               enddo
 
          monthisccphigh(ilon,ilat)=monthisccphigh(ilon,ilat)+               &
                                   isccphighday(ilat,ilon,jour)


          monthisccpliq(ilon,ilat,3)=monthisccpliq(ilon,ilat,3)+             &
                                     isccpliqday(ilat,ilon,jour,3)      
          monthisccpliq(ilon,ilat,4)=monthisccpliq(ilon,ilat,4)+             &
                                     isccpliqday(ilat,ilon,jour,4)

          monthisccpice(ilon,ilat,3)=monthisccpice(ilon,ilat,3)+             &
                                     isccpiceday(ilat,ilon,jour,3)      
          monthisccpice(ilon,ilat,4)=monthisccpice(ilon,ilat,4)+             &
                                     isccpiceday(ilat,ilon,jour,4)
      
do icat=1,catmax
          monthisccpun(ilon,ilat,3,icat)=monthisccpun(ilon,ilat,3,icat)+             &
                                     isccpunday(ilat,ilon,jour,3,icat)      
          monthisccpun(ilon,ilat,4,icat)=monthisccpun(ilon,ilat,4,icat)+             &
                                     isccpunday(ilat,ilon,jour,4,icat)
enddo

          monthcolcloud(ilon,ilat)=monthcolcloud(ilon,ilat)+                 &
                                   colcloudday(ilat,ilon,jour)
       
          isccpdaypermonth(ilon,ilat)=isccpdaypermonth(ilon,ilat)+1
          endif

          if ( isccpinddaylow(ilat,ilon,jour).gt.0 ) then
             do ihist=1,histmax-1
                if( (isccplowday(ilat,ilon,jour).ge.histmod(ihist)) .and.   &
                    (isccplowday(ilat,ilon,jour).lt.histmod(ihist+1)) )then          
                    hlow(ilon,ilat,ihist)=hlow(ilon,ilat,ihist)+1
                endif
             enddo

          monthisccpliq(ilon,ilat,1)=monthisccpliq(ilon,ilat,1)+             &
                                     isccpliqday(ilat,ilon,jour,1)

          monthisccpice(ilon,ilat,1)=monthisccpice(ilon,ilat,1)+             &
                                     isccpiceday(ilat,ilon,jour,1)      
do icat=1,catmax
          monthisccpun(ilon,ilat,1,icat)=monthisccpun(ilon,ilat,1,icat)+             &
                                     isccpunday(ilat,ilon,jour,1,icat)
enddo

          monthisccplow(ilon,ilat)=monthisccplow(ilon,ilat)+                 &
                                   isccplowday(ilat,ilon,jour)
          isccpdaypermonthlow(ilon,ilat)=isccpdaypermonthlow(ilon,ilat)+ 1 
          endif

          if ( isccpinddaymid(ilat,ilon,jour).gt.0 ) then
             do ihist=1,histmax-1
                if( (isccpmidday(ilat,ilon,jour).ge.histmod(ihist)) .and.    &
                    (isccpmidday(ilat,ilon,jour).lt.histmod(ihist+1)) )then           
                    hmid(ilon,ilat,ihist)=hmid(ilon,ilat,ihist)+1
                endif
             enddo

          monthisccpliq(ilon,ilat,2)=monthisccpliq(ilon,ilat,2)+             &
                                     isccpliqday(ilat,ilon,jour,2)

          monthisccpice(ilon,ilat,2)=monthisccpice(ilon,ilat,2)+             &
                                     isccpiceday(ilat,ilon,jour,2)   

do icat=1,catmax
          monthisccpun(ilon,ilat,2,icat)=monthisccpun(ilon,ilat,2,icat)+             &
                                     isccpunday(ilat,ilon,jour,2,icat)
enddo
   
          monthisccpmid(ilon,ilat)=monthisccpmid(ilon,ilat)+                 &
                                   isccpmidday(ilat,ilon,jour)
          isccpdaypermonthmid(ilon,ilat)=isccpdaypermonthmid(ilon,ilat)+ 1 
          endif

 
do nphase=1,4
   isccptemp=isccpiceday(ilat,ilon,jour,nphase)+isccpliqday(ilat,ilon,jour,nphase)
   if(inddayphase(ilat,ilon,jour,nphase).gt.0.)then
      monthisccpphase(ilon,ilat,nphase)=monthisccpphase(ilon,ilat,nphase)+   &
                                        isccpphaseday(ilat,ilon,jour,nphase)
      indmonthphase(ilon,ilat,nphase)=indmonthphase(ilon,ilat,nphase)+1
   else
   isccptemp=sum(isccpunday(ilat,ilon,jour,nphase,:))
      if(isccptemp.gt.0.)then
         indmonthphase2(ilon,ilat,nphase)=indmonthphase2(ilon,ilat,nphase)+1   
      endif
   endif
enddo


        enddo
     enddo
 enddo


!******** CALCULATION OF MONTHLY DIAGNOSTIC WITH MATCHING INDEXES ***********!

 do ilat=1,latmax-1
    do ilon=1,lonmax-1


        if ( indmonthcoltemp(ilon,ilat).ne.0 ) then
           monthcoltemp(ilon,ilat)=                                        &
           monthcoltemp(ilon,ilat)/indmonthcoltemp(ilon,ilat)
        else
           monthcoltemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hcoltemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif

        if ( indmonthhightemp(ilon,ilat).ne.0 ) then
           monthhightemp(ilon,ilat)=                                        &
           monthhightemp(ilon,ilat)/indmonthhightemp(ilon,ilat)
        else
           monthhightemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hhightemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif

        if ( indmonthmidtemp(ilon,ilat).ne.0 ) then
           monthmidtemp(ilon,ilat)=                                        &
           monthmidtemp(ilon,ilat)/indmonthmidtemp(ilon,ilat)
        else
           monthmidtemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hmidtemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif

        if ( indmonthlowtemp(ilon,ilat).ne.0 ) then
           monthlowtemp(ilon,ilat)=                                        &
           monthlowtemp(ilon,ilat)/indmonthlowtemp(ilon,ilat)
        else
           monthlowtemp(ilon,ilat)=-999
           do ihisttemp=1,histmax-1
              hlowtemp(ilon,ilat,ihisttemp)=-999
           enddo
        endif


        if ( isccpdaypermonth(ilon,ilat).ne.0 ) then
           monthisccphigh(ilon,ilat)=                                        &
           monthisccphigh(ilon,ilat)/isccpdaypermonth(ilon,ilat)

          monthisccpliq(ilon,ilat,3)=                                        &
          monthisccpliq(ilon,ilat,3)/isccpdaypermonth(ilon,ilat)
          monthisccpice(ilon,ilat,3)=                                        &
          monthisccpice(ilon,ilat,3)/isccpdaypermonth(ilon,ilat)

          monthisccpliq(ilon,ilat,4)=                                        &
          monthisccpliq(ilon,ilat,4)/isccpdaypermonth(ilon,ilat)
          monthisccpice(ilon,ilat,4)=                                        &
          monthisccpice(ilon,ilat,4)/isccpdaypermonth(ilon,ilat)


           monthcolcloud(ilon,ilat)=                                         &
           monthcolcloud(ilon,ilat)/isccpdaypermonth(ilon,ilat)
           monthcolclear(ilon,ilat)=                                         &
           1-monthcolcloud(ilon,ilat)

do icat=1,catmax
           monthisccpun(ilon,ilat,3,icat)=                                        &
           monthisccpun(ilon,ilat,3,icat)/isccpdaypermonth(ilon,ilat)

           monthisccpun(ilon,ilat,4,icat)=                                        &
           monthisccpun(ilon,ilat,4,icat)/isccpdaypermonth(ilon,ilat)
enddo

           if (indmonthheight(ilon,ilat).ne.0) then
              monthheight(ilon,ilat)=                                        &
              monthheight(ilon,ilat)/indmonthheight(ilon,ilat)
              monthheight2(ilon,ilat)=                                        &
              monthheight2(ilon,ilat)/indmonthheight(ilon,ilat)
           else
              monthheight(ilon,ilat)=-9999.
              monthheight2(ilon,ilat)=-9999.          
           endif
        else
      
           monthisccpliq(ilon,ilat,4)=-9999.
           monthisccpice(ilon,ilat,4)=-9999.
           monthisccpun(ilon,ilat,4,:)=-9999.


           monthisccpliq(ilon,ilat,3)=-9999.
           monthisccpice(ilon,ilat,3)=-9999.
           monthisccpun(ilon,ilat,3,:)=-9999.

           monthisccphigh(ilon,ilat)=-9999.
           monthcolcloud(ilon,ilat)=-9999.
           monthcolclear(ilon,ilat)=-9999.
           monthheight(ilon,ilat)=-9999.
           monthheight2(ilon,ilat)=-9999.          
 
           indtotmean(ilon,ilat)=-999. 
          do ihist=1,histmax2-1
              hheight(ilon,ilat,ihist)=-999.
           enddo
           do ihist=1,histmax-1
              hhigh(ilon,ilat,ihist)=-999.
           enddo
        endif



        if ( isccpdaypermonthlow(ilon,ilat).gt.0 ) then
           monthisccplow(ilon,ilat)=                                         &
           monthisccplow(ilon,ilat)/isccpdaypermonthlow(ilon,ilat)

          monthisccpliq(ilon,ilat,1)=                                        &
          monthisccpliq(ilon,ilat,1)/isccpdaypermonthlow(ilon,ilat)
          monthisccpice(ilon,ilat,1)=                                        &
          monthisccpice(ilon,ilat,1)/isccpdaypermonthlow(ilon,ilat)

do icat=1,catmax
           monthisccpun(ilon,ilat,1,icat)=                                        &
           monthisccpun(ilon,ilat,1,icat)/isccpdaypermonthlow(ilon,ilat)
enddo

        else
           monthisccpice(ilon,ilat,1)=-9999.
           monthisccpliq(ilon,ilat,1)=-9999.
           monthisccpun(ilon,ilat,1,:)=-9999.
           monthisccplow(ilon,ilat)=-9999.

           do ihist=1,histmax-1
              hlow(ilon,ilat,ihist)=-999.
           enddo
        endif


        if ( isccpdaypermonthmid(ilon,ilat).gt.0 ) then
           monthisccpmid(ilon,ilat)=                                         &
           monthisccpmid(ilon,ilat)/isccpdaypermonthmid(ilon,ilat)

          monthisccpliq(ilon,ilat,2)=                                        &
          monthisccpliq(ilon,ilat,2)/isccpdaypermonthmid(ilon,ilat)
          monthisccpice(ilon,ilat,2)=                                        &
          monthisccpice(ilon,ilat,2)/isccpdaypermonthmid(ilon,ilat)

do icat=1,catmax
           monthisccpun(ilon,ilat,2,icat)=                                        &
           monthisccpun(ilon,ilat,2,icat)/isccpdaypermonthmid(ilon,ilat)
enddo  

        else
           monthisccpice(ilon,ilat,2)=-9999.
           monthisccpliq(ilon,ilat,2)=-9999.
           monthisccpun(ilon,ilat,2,:)=-9999.
           monthisccpmid(ilon,ilat)=-9999.

           do ihist=1,histmax-1
              hmid(ilon,ilat,ihist)=-999.
           enddo

        endif

        do nphase=1,4
           if ( indmonthphase(ilon,ilat,nphase).gt.0)then
              monthisccpphase(ilon,ilat,nphase)=                             &
              monthisccpphase(ilon,ilat,nphase)/indmonthphase(ilon,ilat,nphase)
           else
              if (indmonthphase2(ilon,ilat,nphase).gt.0)then
                 monthisccpphase(ilon,ilat,nphase)=-777.  
              else
                 monthisccpphase(ilon,ilat,nphase)=-9999.
              endif
           endif
        enddo

                                    
    enddo
enddo



!***************************** SAVE THE MAP FILES ***************************!

! MapLowMidHigh_

file8=trim(file6)//'.nc'    ! name of output ncdf map file
file9=trim(file3(25:55))    ! period of map file (description of ncdf file)

call create_mapnc(file8,file9,lonmid,latmid,resd,dimidsm,gcm,lonmax-1,latmax-1)
call map_recvar2nc2(monthisccplow,monthisccpmid,monthisccphigh,monthcolcloud,&
                     monthcolclear,dimidsm,file8,lonmax-1,latmax-1)

! MapHigh_

file8=trim(file66)//'.nc'    ! name of output ncdf map file

call create_maphighnc(file8,file9,lonmid,latmid,resd,dimidsm,gcm,lonmax-1,latmax-1)
call maphigh(monthisccphigh,monthheight,monthheight2,dimidsm,file8,lonmax-1,latmax-1)

! Change NaN value from -9999 to -999 to fit with the GEWEX standard
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthisccplow(ilon,ilat)==-9999.)
monthisccplow(ilon,ilat)=-999
endforall
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthisccpmid(ilon,ilat)==-9999.)
monthisccpmid(ilon,ilat)=-999
endforall
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthisccphigh(ilon,ilat)==-9999.)
monthisccphigh(ilon,ilat)=-999
endforall
forall(ilon=1:lonmax-1, ilat=1:latmax-1, monthcolcloud(ilon,ilat)==-9999.)
monthcolcloud(ilon,ilat)=-999
endforall


! MapLowMidHigh_Phase

file8=trim(file11)//'.nc'    ! name of output ncdf map file
call create_mapnc_phase(file8,file9,lonmid,latmid,resd,dimidsm,dimidsm2,gcm,lonmax-1,latmax-1)
call map_recvar2nc2phaseocc2(monthisccpliq,monthisccpice,monthisccpun,        &
                              monthisccpphase,dimidsm,dimidsm2,file8,          &
                              lonmax-1,latmax-1, catmax)


! Deallocate daily & monthly map variables
print *, 'deallocate daily & monthly map variables'

! fixme : here we are in the "case("lmdz")" section of the select(model). 
! model is always == 'lmdz' (sauf erreur !)
if(model=='lmdz')then
  deallocate(hlow,hmid,hhigh,hcol,hheight)
  deallocate(monthisccplow,monthisccpmid,monthisccphigh)
  deallocate(monthcolcloud,isccpdaypermonth)
  deallocate(monthisccpliq,monthisccpice,monthisccpun)
  deallocate(indmonthphase,monthisccpphase)
  deallocate(indmonthphase2)
  deallocate(inddayphase)
  deallocate(isccplowday,isccpmidday,isccphighday)
  deallocate(isccpliqday,isccpiceday,isccpunday)
  deallocate(colcloudday,isccpindday)
  deallocate(isccpinddaylow,isccpinddaymid)
  deallocate(isccpdaypermonthlow,isccpdaypermonthmid)
  deallocate(indtotmean,indtot)
  deallocate(indmonthheight,monthheight,indheight,heightday)
  deallocate(monthheight2,heightday2)
  deallocate(lowtemp,midtemp,hightemp,coltemp)
  deallocate(indlowtemp,indmidtemp,indhightemp,indcoltemp)
  deallocate(hlowtemp,hmidtemp,hhightemp,hcoltemp)

endif

print *, 'map file recorded'




!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!! PART III : CLOUDY MAP3D FILES !!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


! Allocation / initialization of MAP3D monthly variables
print *, 'allocation / initialization of MAP3D monthly variables'
   allocate(indphasepermonth(lonmax-1,latmax-1,altmax))

   allocate(indpermonth(lonmax-1,latmax-1,altmax))

   allocate(monthcloudfract(lonmax-1,latmax-1,altmax),                           &
            monthclearfract(lonmax-1,latmax-1,altmax)) 

   allocate(monthuncertfract(lonmax-1,latmax-1,altmax))


   allocate(monthicecloud(lonmax-1,latmax-1,altmax),                       &
            monthwatercloud(lonmax-1,latmax-1,altmax),                       &
            indphasemonth(lonmax-1,latmax-1,altmax))
   allocate(monthuncloud(lonmax-1,latmax-1,altmax,catmax))
   allocate(monthphasecloud(lonmax-1,latmax-1,altmax))
   allocate(indmonthphase3D(lonmax-1,latmax-1,altmax))
   allocate(indphasefractday(latmax-1,lonmax-1,altmax,daymax))
   allocate(monthcftemp(lonmax-1,latmax-1,tempmax-1))
   allocate(monthcftempice(lonmax-1,latmax-1,tempmax-1))
   allocate(monthcftempliq(lonmax-1,latmax-1,tempmax-1))
   allocate(monthcftempphase(lonmax-1,latmax-1,tempmax-1))
   allocate(indmonthphasetemp(lonmax-1,latmax-1,tempmax-1))
   allocate(indcftemppermonth(lonmax-1,latmax-1,tempmax-1))
allocate(cftempphaseday(latmax-1,lonmax-1,tempmax-1,daymax))


 cftempphaseday(:,:,:,:)=0;
indpermonth(:,:,:)=0;indphasepermonth(:,:,:)=0;
indmonthphase3D(:,:,:)=0; monthphasecloud(:,:,:)=0;
monthcloudfract(:,:,:)=0;monthclearfract(:,:,:)=0;
monthuncertfract(:,:,:)=0;
indphasemonth(:,:,:)=0;monthicecloud(:,:,:)=0;monthwatercloud(:,:,:)=0;
monthuncloud(:,:,:,:)=0;
indphasefractday(:,:,:,:)=0;
monthcftemp(:,:,:)=0
monthcftempice(:,:,:)=0
monthcftempliq(:,:,:)=0
monthcftempphase(:,:,:)=0
indmonthphasetemp(:,:,:)=0
indcftemppermonth(:,:,:)=0


!! OCCURRENCES FILE IN 3D WITH TEMPERATURE
  do jour=1,31
    do itemp=1,tempmax-1
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indcftemp(ilat,ilon,itemp,jour).gt.0) then
     indcftemppermonth(ilon,ilat,itemp) = indcftemppermonth(ilon,ilat,itemp)+indcftemp(ilat,ilon,itemp,jour)
     monthcftemp(ilon,ilat,itemp)=monthcftemp(ilon,ilat,itemp)   &
                                +cftempday(ilat,ilon,itemp,jour)
     monthcftempice(ilon,ilat,itemp)=monthcftempice(ilon,ilat,itemp)   &
                                +cftempiceday(ilat,ilon,itemp,jour)
     monthcftempliq(ilon,ilat,itemp)=monthcftempliq(ilon,ilat,itemp)   &
                                +cftempliqday(ilat,ilon,itemp,jour)
          endif
         enddo
       enddo
    enddo
  enddo


file8=trim(file13)//'_occ.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

print *, 'MAP3D PHASE files recorded'


 call create_temp3d(file8,file9,lonmid,latmid,tempmid,tempmod_bound,resd,       &
                    dimidsp,tempmax-1,lonmax-1,latmax-1)


print *, 'MAP3D create temp'

 call temp_recvar2nc(monthcftemp,monthcftempliq,monthcftempice,        &
                     indcftemppermonth,dimidsp,file8,tempmax-1,lonmax-1,latmax-1)

monthcftemp(:,:,:)=0
monthcftempice(:,:,:)=0
monthcftempliq(:,:,:)=0
indcftemppermonth(:,:,:)=0


!********* CALCULATION OF DAILY DIAGNOSTIC WITH MATCHING INDEXES ************!

 cfsumtemp=0.

  do jour=1,31
    do itemp=1,tempmax-1
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indcftemp(ilat,ilon,itemp,jour).gt.0) then


 cfsumtemp=cftempliqday(ilat,ilon,itemp,jour)+  &
          cftempiceday(ilat,ilon,itemp,jour)

if(cfsumtemp.gt.0)then
   if(cftempiceday(ilat,ilon,itemp,jour).eq.0.)then
      cftempphaseday(ilat,ilon,itemp,jour)=0.
      indcftempphase(ilat,ilon,itemp,jour)=  &
      indcftempphase(ilat,ilon,itemp,jour)+1
   else
      cftempphaseday(ilat,ilon,itemp,jour)=    &
      cftempiceday(ilat,ilon,itemp,jour)/cfsumtemp
      indcftempphase(ilat,ilon,itemp,jour)=  &
      indcftempphase(ilat,ilon,itemp,jour)+1
   endif
endif

 cftempday(ilat,ilon,itemp,jour)=   &
 cftempday(ilat,ilon,itemp,jour)/indcftemp(ilat,ilon,itemp,jour)
 cftempliqday(ilat,ilon,itemp,jour)=   &
 cftempliqday(ilat,ilon,itemp,jour)/indcftemp(ilat,ilon,itemp,jour)
 cftempiceday(ilat,ilon,itemp,jour)=   &
 cftempiceday(ilat,ilon,itemp,jour)/indcftemp(ilat,ilon,itemp,jour)

if(cftempday(ilat,ilon,itemp,jour).gt.1.)then
print *, ilat,ilon,itemp,jour
print *, cftempday(ilat,ilon,itemp,jour)
print *, cftempliqday(ilat,ilon,itemp,jour)
print *, cftempiceday(ilat,ilon,itemp,jour)

endif

           endif
         enddo
       enddo
    enddo
  enddo


allocate(tot_ind(lonmax-1,latmax-1,altmax), &
         cloud_ind(lonmax-1,latmax-1,altmax), &
         ice_ind(lonmax-1,latmax-1,altmax), &
         water_ind(lonmax-1,latmax-1,altmax), &
         un_ind(lonmax-1,latmax-1,altmax,catmax))


print *, 'test'

    do ialt=1,altmax
       do ilon=1,lonmax-1
         do ilat=1,latmax-1

    cloud_ind(ilon,ilat,ialt)=sum(cloudfractday(ilat,ilon,ialt,:))
    tot_ind(ilon,ilat,ialt)=sum(indday(ilat,ilon,ialt,:))
    ice_ind(ilon,ilat,ialt)=sum(icecloudfractday(ilat,ilon,ialt,:))
    water_ind(ilon,ilat,ialt)=sum(watercloudfractday(ilat,ilon,ialt,:))

         do icat=1,catmax
    un_ind(ilon,ilat,ialt,icat)=sum(uncloudfractday(ilat,ilon,ialt,:,icat))
         enddo

       enddo
    enddo
  enddo



file8=trim(file10)//'_occ.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

print *, 'titi'
 call create_ind3d(file8,file9,lonmid,latmid,altmid,altmod_bound,resd,   &
                       dimidsp,dimidsp2,altmax,lonmax-1,latmax-1)

print *, 'titi2'

 call record_ind3d(cloud_ind,tot_ind,ice_ind,water_ind,un_ind, &
                  dimidsp,dimidsp2,file8,altmax,lonmax-1,latmax-1)
print *, 'titi3'

deallocate(tot_ind,cloud_ind,ice_ind,water_ind,un_ind)


isccptemp=0.
  do jour=1,31
    do ialt=1,altmax
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indday(ilat,ilon,ialt,jour).gt.0) then
! Cloud fraction / Clear fraction / Saturated fraction, daily
         cloudfractday(ilat,ilon,ialt,jour) =                                &
         cloudfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)        
         clearfractday(ilat,ilon,ialt,jour) =                                &
         clearfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
         uncertfractday(ilat,ilon,ialt,jour) =                               &
         uncertfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
         isccptemp=icecloudfractday(ilat,ilon,ialt,jour)+                     &
                   watercloudfractday(ilat,ilon,ialt,jour)
if(isccptemp.gt.0)then
   if(icecloudfractday(ilat,ilon,ialt,jour).eq.0)then
   phasefractday(ilat,ilon,ialt,jour) = 0.
   indphasefractday(ilat,ilon,ialt,jour) = indphasefractday(ilat,ilon,ialt,jour)+1
   else
   indphasefractday(ilat,ilon,ialt,jour) = indphasefractday(ilat,ilon,ialt,jour)+1
   phasefractday(ilat,ilon,ialt,jour) =                                        &
   icecloudfractday(ilat,ilon,ialt,jour)/isccptemp
   endif
else
         isccptemp=sum(uncloudfractday(ilat,ilon,ialt,jour,:))
         if(isccptemp.gt.0)then
           phasefractday(ilat,ilon,ialt,jour)=-777. 
         endif
endif

         icecloudfractday(ilat,ilon,ialt,jour) =                                &
         icecloudfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)
         watercloudfractday(ilat,ilon,ialt,jour) =                                &
         watercloudfractday(ilat,ilon,ialt,jour)/indday(ilat,ilon,ialt,jour)

do icat=1,catmax
         uncloudfractday(ilat,ilon,ialt,jour,icat) =                                &
         uncloudfractday(ilat,ilon,ialt,jour,icat)/indday(ilat,ilon,ialt,jour)
enddo


          endif

         enddo
       enddo
    enddo
  enddo



!************ INCREMENTING THE MONTHLY DIAGNOSTIC AND INDEXES **************!


  do jour=1,31
    do itemp=1,tempmax-1
       do ilon=1,lonmax-1
         do ilat=1,latmax-1
           if (indcftemp(ilat,ilon,itemp,jour).gt.0) then
     indcftemppermonth(ilon,ilat,itemp) = indcftemppermonth(ilon,ilat,itemp)+1
     monthcftemp(ilon,ilat,itemp)=monthcftemp(ilon,ilat,itemp)   &
                                +cftempday(ilat,ilon,itemp,jour)
     monthcftempice(ilon,ilat,itemp)=monthcftempice(ilon,ilat,itemp)   &
                                +cftempiceday(ilat,ilon,itemp,jour)
     monthcftempliq(ilon,ilat,itemp)=monthcftempliq(ilon,ilat,itemp)   &
                                +cftempliqday(ilat,ilon,itemp,jour)

!indcftempphase
if(indcftempphase(ilat,ilon,itemp,jour).gt.0.)then
    monthcftempphase(ilon,ilat,itemp) = monthcftempphase(ilon,ilat,itemp)+  &
                                      cftempphaseday(ilat,ilon,itemp,jour)
    indmonthphasetemp(ilon,ilat,itemp) = indmonthphasetemp(ilon,ilat,itemp)+1 
endif

           endif
         enddo
       enddo
    enddo
  enddo

  do jour=1,31
    do ialt=1,altmax
      do ilat=1,latmax-1
       do ilon=1,lonmax-1



  if (indday(ilat,ilon,ialt,jour).gt.0) then
     indpermonth(ilon,ilat,ialt) = indpermonth(ilon,ilat,ialt)+1
    indphasepermonth(ilon,ilat,ialt) = indphasepermonth(ilon,ilat,ialt)+ &
                                      indday(ilat,ilon,ialt,jour)

     ! Monthly Fraction 
     monthcloudfract(ilon,ilat,ialt) = monthcloudfract(ilon,ilat,ialt) +     &
                                       cloudfractday(ilat,ilon,ialt,jour)
     monthclearfract(ilon,ilat,ialt) = monthclearfract(ilon,ilat,ialt) +     &
                                       clearfractday(ilat,ilon,ialt,jour)
  
     monthuncertfract(ilon,ilat,ialt) = monthuncertfract(ilon,ilat,ialt) +   &
                                        uncertfractday(ilat,ilon,ialt,jour)

    monthicecloud(ilon,ilat,ialt) = monthicecloud(ilon,ilat,ialt) +     &
                                    icecloudfractday(ilat,ilon,ialt,jour)
    monthwatercloud(ilon,ilat,ialt) = monthwatercloud(ilon,ilat,ialt) +     &
                                    watercloudfractday(ilat,ilon,ialt,jour)

do icat=1,catmax
    monthuncloud(ilon,ilat,ialt,icat) = monthuncloud(ilon,ilat,ialt,icat) +     &
                                    uncloudfractday(ilat,ilon,ialt,jour,icat)
enddo

!!! changement ici 
!! a terminer la phase 3D
if(indphasefractday(ilat,ilon,ialt,jour).gt.0.)then
    monthphasecloud(ilon,ilat,ialt) = monthphasecloud(ilon,ilat,ialt) +    &
                                      phasefractday(ilat,ilon,ialt,jour)
    indmonthphase3D(ilon,ilat,ialt) = indmonthphase3D(ilon,ilat,ialt)+1 
endif

  endif



       enddo
      enddo
    enddo
  enddo


    do ialt=1,altmax
      do ilat=1,latmax-1
       do ilon=1,lonmax-1

  if (indphasepermonth(ilon,ilat,ialt).lt.sum(monthuncloud(ilon,ilat,ialt,:))) then

print *, ilon,ilat,ialt,indphasepermonth(ilon,ilat,ialt),monthuncloud(ilon,ilat,ialt,:)
  endif

enddo
enddo
enddo


! Deallocate daily MAP3D variables
print *, 'deallocate daily MAP3D variables'

if(model=='lmdz')then
 deallocate(indday)
 deallocate(cloudfractday, clearfractday,uncertfractday)
 deallocate(indphaseday,icecloudfractday,watercloudfractday)
 deallocate(uncloudfractday)
 deallocate(cftempday,cftempliqday,cftempiceday,indcftemp,indcftempphase)

endif


!******** CALCULATION OF MONTHLY DIAGNOSTIC WITH MATCHING INDEXES ***********!

do ilat=1,latmax-1
    do ilon=1,lonmax-1

       do itemp=1,tempmax-1
          if(indmonthphasetemp(ilon,ilat,itemp).gt.0)then
          monthcftempphase(ilon,ilat,itemp)=                                 &
          monthcftempphase(ilon,ilat,itemp)/indmonthphasetemp(ilon,ilat,itemp)
          else
          monthcftempphase(ilon,ilat,itemp)=-9999.   
          endif

          if(indcftemppermonth(ilon,ilat,itemp).gt.0 )then
          monthcftemp(ilon,ilat,itemp)=    &
          monthcftemp(ilon,ilat,itemp)/indcftemppermonth(ilon,ilat,itemp)
          monthcftempice(ilon,ilat,itemp)=    &
          monthcftempice(ilon,ilat,itemp)/indcftemppermonth(ilon,ilat,itemp)
          monthcftempliq(ilon,ilat,itemp)=    &
          monthcftempliq(ilon,ilat,itemp)/indcftemppermonth(ilon,ilat,itemp)

          else
          monthcftemp(ilon,ilat,itemp)=-9999.   
          monthcftempliq(ilon,ilat,itemp)=-9999.   
          monthcftempice(ilon,ilat,itemp)=-9999.   

          endif
       enddo


       do ialt=1,altmax    


        if ( indmonthphase3D(ilon,ilat,ialt).gt.0 ) then
        monthphasecloud(ilon,ilat,ialt)=                                    &
        monthphasecloud(ilon,ilat,ialt)/indmonthphase3D(ilon,ilat,ialt)  
        else
        monthphasecloud(ilon,ilat,ialt)=-9999.     
           do jour = 1,31
              if(phasefractday(ilat,ilon,ialt,jour).eq.-777.)then
              monthphasecloud(ilon,ilat,ialt)=-777.
              endif
           enddo
           
        endif

        if ( indpermonth(ilon,ilat,ialt).gt.0 ) then
        monthcloudfract(ilon,ilat,ialt)=                                     &
        monthcloudfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
        monthclearfract(ilon,ilat,ialt)=                                     &
        monthclearfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
        monthuncertfract(ilon,ilat,ialt)=                                    &
        monthuncertfract(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)

        monthicecloud(ilon,ilat,ialt)=                                    &
        monthicecloud(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
        monthwatercloud(ilon,ilat,ialt)=                                    &
        monthwatercloud(ilon,ilat,ialt)/indpermonth(ilon,ilat,ialt)
do icat=1,catmax
        monthuncloud(ilon,ilat,ialt,icat)=                                    &
        monthuncloud(ilon,ilat,ialt,icat)/indpermonth(ilon,ilat,ialt)
enddo

        else
        monthcloudfract(ilon,ilat,ialt)=-9999.
        monthclearfract(ilon,ilat,ialt)=-9999.
      
        monthuncertfract(ilon,ilat,ialt)=-9999.
        monthuncloud(ilon,ilat,ialt,:)=-9999.
        monthwatercloud(ilon,ilat,ialt)=-9999.
        monthicecloud(ilon,ilat,ialt)=-9999. 
        indphasepermonth(ilon,ilat,ialt)=-9999. 
       endif


      enddo
    enddo
  enddo


deallocate(phasefractday,indphasefractday)
! Deallocate monthly indexes

!**************************** SAVE THE MAP3D FILES **************************!

! name of output netcdf MAP3D file
file8='Map3D330m_'//trim(file5)//'_'//trim(switch)//'_'//trim(gcm)//'_'//trim(switch2)//'_'//trim(version)//'.nc'
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 



 call create_profnc(file8,file9,lonmid,latmid,altmid,altmod_bound,resd,       &
                    dimidsp,altmax,lonmax-1,latmax-1)


 call prof_recvar2nc(monthcloudfract,monthclearfract,monthuncertfract,        &
                     dimidsp,file8,altmax,lonmax-1,latmax-1)

print *, 'MAP3D files recorded'

file8=trim(file10)//'.nc'   ! name of output netcdf MAP3D file FIXME
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

 call create_depolnc3d(file8,file9,lonmid,latmid,altmid,altmod_bound,resd,   &
                       dimidsp,dimidsp2,altmax,lonmax-1,latmax-1)


 call depol_recvar2ncocc(monthicecloud,monthwatercloud,monthuncloud,          &
                         monthphasecloud,indphasepermonth,dimidsp,dimidsp2,   &
                         file8,altmax,lonmax-1,latmax-1)


! check the allocation 
 if (OK_buffer/=0) print *,'--- buffer allocation error '   

  deallocate(monthphasecloud,stat = OK_buffer)
if (OK_buffer/=0) print *,'--- buffer allocation error ' 
  deallocate(indphasepermonth,stat = OK_buffer)
if (OK_buffer/=0) print *,'--- buffer allocation error ' 


file8=trim(file13)//'.nc'   ! name of output netcdf MAP3D file
file9=trim(file3(25:55))   ! period of MAP3D file (description of ncdf file) 

print *, 'MAP3D PHASE files recorded'


 call create_temp3d(file8,file9,lonmid,latmid,tempmid,tempmod_bound,resd,       &
                    dimidsp,tempmax-1,lonmax-1,latmax-1)


print *, 'MAP3D create temp'

 call temp_recvar2nc(monthcftemp,monthcftempliq,monthcftempice,        &
                     monthcftempphase,dimidsp,file8,tempmax-1,lonmax-1,latmax-1)

 
print *, 'MAP3D record temp'


! Deallocate daily & monthly diagSR variables
print *, 'deallocate daily & monthly MAP3D variables'

if(model=='lmdz')then
  deallocate(monthcloudfract,monthclearfract,monthuncertfract)
  deallocate(monthicecloud,monthwatercloud,indphasemonth,indpermonth)
  deallocate(monthuncloud)
  deallocate(indcftemppermonth)
deallocate(indmonthphasetemp)
deallocate(monthcftempphase)

  deallocate(monthcftemp)
deallocate(monthcftempliq)
deallocate(monthcftempice)

endif






!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!! PART II : DIAGSR FILES !!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!


! Allocation / initialization of diagSR monthly variables
print *, 'allocation / initialization of diagSR monthly variables'

file8=trim(file7)//'.nc'   ! name of output netcdf diagSR file
file9=trim(file3(25:55))   ! period of diagSR file (description of ncdf file)

print *, 'allocation monthdiagSR terminé'


print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,1),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,2),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,3),4),3),2),1)

forall(ilon=1:lonmax-1, ilat=1:latmax-1,  ialt=1:altmax, indnan(ilat,ilon,ialt) == 0)
 diagSR(ilon,ilat,ialt,1:diagmax-1)=-9999.
 diagSRpha(ilon,ilat,ialt,1:diagmax-8,:)=-9999.

endforall



 call create_diagnc(file8,file9,lonmid,latmid,altmid,altmod_bound,srmod,resd,dimidsd,dimidsdb, &
                    altmax,lonmax-1,latmax-1)

print *, 'creation fichier diag final'

 call diag_recvar2nc3(diagSR,dimidsd,dimidsdb,file8,altmax,lonmax-1,latmax-1)


file8=trim(file12)//'.nc'   ! name of output netcdf diagSR file
file9=trim(file3(25:55))   ! period of diagSR file (description of ncdf file)

print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,1),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,2),4),3),2),1)
print *, sum(sum(sum(sum(diagSRpha(:,:,:,:,3),4),3),2),1)

 call create_diagncpha(file8,file9,lonmid,latmid,altmid,altmod_bound,srmod,resd,dimidsd,dimidsdb, &
                    altmax,lonmax-1,latmax-1)

print *, 'creation fichier diag final'

 call diag_recvar2nc3pha(diagSRpha,dimidsd,dimidsdb,file8,altmax,lonmax-1,latmax-1)






! Deallocate daily & monthly diagSR variables
print *, 'deallocate daily & monthly diagSR variables'

if(model=='lmdz')then
   deallocate(diagSR,diagSRpha)!,monthdiagSR15,monthdiagSR1)
endif
! Deallocate daily & monthly diagSR variables
print *, 'deallocate daily & monthly diagSR variables'



621 continue



if(model=='lmdz')then
  deallocate(latmod,lonmod,prestop,altmod,srmod,pr2mod,atbrmod,srdepmod,depolmod,lonmid,latmid,altmid,tempmod,stat = OK_buffer) !crmod
endif

 close(1)

CASE ("chimere")
  continue
 
CASE ("wrf")
  deallocate(latmod,lonmod,prestop,altmod,srmod,lonmid,latmid,altmid,stat = OK_buffer) !crmod
  close(1)

CASE DEFAULT
  print *, "error" 

END SELECT

 102  format(f10.2,3e13.5)
 103  format((2x,f10.2),5(2x,e13.5))
 104 format(2(2x,I5),(2x,e13.5),4(2x,f8.2))
 106  format(3(f8.3,1x),2(e12.5,1x),7(f8.3,1x))
 107  format(9(f8.2,1x))
 108 format(4(f8.3,1x),(1x,e13.6))

print *, 'END OF PROGRAM'
print *, '-----------------------------------------------------------------------'
print *, '-----------------------------------------------------------------------'






contains


!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!

!****************************************************************************!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!*                                                                          *!
!*                            END OF PROGRAM                                *!
!*                                                                          *!
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!
!****************************************************************************!

!****************************************************************************!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************!





!****************************************************************************!
!****************************************************************************!
!**************************    SUBROUTINE LIST    ***************************!
!****************************************************************************!
!****************************************************************************!

! reading routines put the data directly into variables defined in the main program.
! extracting them into their own file will require a common block (urg)

!----------------------------------------------------------------------------!
! **** SDSREAD8 ****  This routine allows to read 8-bit SDS-variable on hdf  !
!                     files                                                  !
!----------------------------------------------------------------------------!
! filename  : calipso file name
! varname   : SDS 8-bit name variable to extract from calipso file
!----------------------------------------------------------------------------!
! var       : SDS 8-bit variable to extract from calipso file
!----------------------------------------------------------------------------!
! i         : loop on the indice of calipso variables
! name      : SDS 8-bit name variable to extract from calipso file
! ret       : 
! start     :
! stride    :
! edges     :
! sd_id     : some id
! sds_id    : some id
! OK_buffer :
! npts      : value of one dimension of the read var
! nprofs    : value of one dimension of the read var
! nb        : some loop index
! index     :
! rank      :
! istat     :
! attributes:
! num_type  :
!                                                                            !
!----------------------------------------------------------------------------!
subroutine sdsread8(var,filename,varname,ret)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: i, ret, start(2) = 0, stride(2) = 1 ,edges(2)
	integer 	:: sd_id, sds_id, OK_buffer
        integer         :: npts, nprofs, nb
        real*8,dimension(:,:),allocatable :: var
	integer*4	:: index, rank, istat, attributes, num_type
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
        character       :: varname*100
    
	sd_id = sfstart (filename,	DFACC_READ)

! loop on the variables number (<100)
do i=0,100

   nb=i
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)
                       
        if(name==trim(varname))exit
          
 enddo

       ! Retrieving the dimesions & allocating memorie
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)

      allocate(var(npts, nprofs),stat = OK_buffer)
      if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(name),'error' 

        edges = [npts, nprofs]
        ret = sfrdata (sds_id, start, stride, edges, var)

	! Do something with the data
      if (ret.eq.-1) then
		print *,'ERROR'
goto 555
      end if

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)
555 continue


end subroutine sdsread8
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! **** SDSREAD ****  This routine allows to read a 2dim SDS-variable on hdf  !
!                    files                                                   !
!----------------------------------------------------------------------------!
subroutine sdsread(var,filename,varname)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: i, ret, start(2) = 0, stride(2) = 1 ,edges(2)
	integer 	:: sd_id, sds_id
        integer         :: npts, nprofs, nb, OK_buffer
        real,dimension(:,:),allocatable :: var
	integer*4	:: index, rank, istat, attributes, num_type
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
         character       :: varname*100
 
    
	sd_id = sfstart (filename,	DFACC_READ)
  

do i=2,100

   nb=i
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)

        if(name==trim(varname))exit

enddo

    
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)
      allocate(var(npts, nprofs),stat = OK_buffer)
      if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(name),'error' 
        edges = [npts, nprofs]
        ret = sfrdata (sds_id, start, stride, edges, var)


	! Do something with the data
      if (ret.eq.-1) then
		print *,'ERROR'
      end if

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)

end subroutine sdsread
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
!*****  METAREAD ***** This routine allows to read a Meta-variable on hdf    !
!                      files                                                 !
!----------------------------------------------------------------------------!
subroutine metaread(var,varname,filename)
 
	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: ret, istat, OK_buffer
	integer		:: file_id, vdata_ref, vdata_id

  character :: varname*30
	character	:: filename*1024

	
 	real*4,dimension(:),allocatable  ::  var	
  
! Allocate var size
   if(varname=='Lidar_Data_Altitudes')then
      allocate(var(583), stat=OK_buffer)
   if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(varname),'error' 
        else
      allocate(var(33), stat=OK_buffer)
   if (OK_buffer/=0) print *,'--- buffer allocation of ',trim(varname),'error' 
        endif



	file_id = hopen (filename,	DFACC_READ, 0)
 
	! initialize vdata interface,find vdata called 'metadata '
	! attach to this vdata
	istat = vfstart (file_id)
	vdata_ref = vsffnd (file_id, 'metadata')
	vdata_id = vsfatch (file_id, vdata_ref, 'r')

	! reads Varname

	istat = vsfsfld (vdata_id, varname)
	ret = vsfread (vdata_id, var, 1, FULL_INTERLACE)
	if (ret.ne.-1) then
           continue
        else

	end if
 
 	! TRES IMPORTANT : il faut appeler vsfseek entre chaque vsfsfld
	! sinon la lecture ne fonctionne pas !
 
	ret = vsfseek (vdata_id, 0)
 
	! Close all the crap
	istat = vsfdtch(vdata_id)
	istat = vfend (file_id)
	istat = hclose (file_id)
 
	if (istat.eq.0) then
		print *,'Reading HDF File OK.'
	else
		print *,'ERROR READING HDF'
        endif

end subroutine metaread
!----------------------------------------------------------------------------!


!----------------------------------------------------------------------------!
! ***** NPROF ***** This routine find the hdf file number of profiles        !
!----------------------------------------------------------------------------!
subroutine nprof(filename,nb,nprofs)

	implicit none
	include "hdf.f90"
	include "dffunc.f90"
 
	integer 	:: ret
	integer 	:: sd_id, sds_id
        integer         :: npts, nprofs, nb
	integer*4	:: rank, istat, attributes, num_type !index
	integer*4 	:: dim_sizes (32)
	character(len=232) :: name, filename
 
     
	sd_id = sfstart (filename,	DFACC_READ)
  
	! Select and read the var
	sds_id = sfselect (sd_id, nb)
        istat = sfginfo (sds_id, name, rank, dim_sizes, num_type, attributes)
        
        ! Rtrieving the dimesions
        npts = dim_sizes(1)
        nprofs = dim_sizes(2)
print *, 'le nombre de profil du fichier est :', nprofs

	! Close the crap
	ret = sfendacc(sds_id)
	ret = sfend(sd_id)

end subroutine nprof
!----------------------------------------------------------------------------!




end program calmdz

!****************************************************************************!
!                                                                            !
!                               END PROGRAM                                  !
!                                                                            !
!****************************************************************************!


