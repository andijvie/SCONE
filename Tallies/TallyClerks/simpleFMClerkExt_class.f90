module simpleFMClerkExt_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables
  use genericProcedures,          only : fatalError
  use display_func,               only : statusMsg
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill
  use tallyResult_class,          only : tallyResult

  ! Nuclear Data
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Response
  use macroResponse_class,        only : macroResponse

  ! TODO: output final FM

  implicit none
  private

  !!
  !! Simple 1-D fission matrix
  !!
  !! This is a prototype implementation
  !! Uses collision estimator only
  !! Contains only a single map for discretisation
  !!
  !! Notes:
  !!    -> If collision particle has invalid nuclear data type collision is ignored
  !!    -> Collisions in non-fissile materials are ignored
  !!    -> FM is stored in column-major order [prodBin, startBin]
  !!
  !! Private Members:
  !!   map      -> Map to divide phase-space into bins
  !!   resp     -> Response for transfer function (nuFission by default)
  !!   startWgt -> Starting Weigths in each bin
  !!   N        -> Number of Bins
  !!
  !! Interface:
  !!   tallyClerk Interface
  !!
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type simpleFMClerk;
  !!      map { <TallyMapDef> }
  !!  }
  !!
  type, public, extends(tallyClerk) :: simpleFMClerkExt
    private
    !! Map defining the discretisation
    class(tallyMap), allocatable :: map
    type(macroResponse)          :: resp
    integer(shortInt)            :: N = 0 ! Number of bins
    integer(shortInt)            :: window ! NEW: No. of cycles to save
    logical(defBool)             :: doDebug ! NEW: extra prints
    logical(defBool)             :: forceOne ! NEW: forces a homogeneous eigenvector

    ! NEW: Fundamental eigenvector for scaling particles
    real(defReal),dimension(:),allocatable       :: eigVec   

    ! NEW: Stores NOT-normalized FM for multiple cycles
    ! This is a queue (FiLo)
    real(defReal), dimension(:,:,:), allocatable :: tallyMatrix

    ! NEW: Stores normalized FM
    real(defReal), dimension(:,:), allocatable   :: matrix

    ! NEW: Stores the cumilative number of neutrons in the starting bin for multiple cycles
    ! This is a queue (FiLo)
    real(defReal),dimension(:,:),allocatable     :: startWgt

    ! Settings
    logical(defBool) :: handleVirtual = .true.

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleStart
    procedure  :: reportInColl
    procedure  :: closeCycle

    ! Overwrite default run-time result procedure
    procedure  :: getResult

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill

    ! NEW: Solve for the FM eigenvector
    procedure  :: solve 

  end type simpleFMClerkExt

  !!
  !! NEW: FM eigenvector result class
  !!   Stored in column first order
  !!    dim1 -> target bin
  !!    dim2 -> orgin bin
  !!
  type,public, extends( tallyResult) :: FMresult
    integer(shortInt)                       :: N  = 0 ! Size of FM
    real(defReal), dimension(:),allocatable :: eigVec ! FM eigenvector
  end type FMResult


  !!
  !!  NAMECHANGE: Fission matrix result class
  !!   Stored in column first order
  !!    dim1 -> target bin
  !!    dim2 -> orgin bin
  !!    dim3 -> 1 is values; 2 is STDs
  !!
  type,public, extends( tallyResult) :: FMmatrix
    integer(shortInt)                            :: N  = 0 ! Size of FM
    real(defReal), dimension(:,:,:), allocatable :: FM     ! FM proper
  end type FMmatrix

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(simpleFMClerkExt), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    character(nameLen), intent(in)      :: name

    ! Assign name
    call self % setName(name)

    ! Read map
    call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Read size of the map
    self % N = self % map % bins(0)

    ! NEW: read size of window
    call dict % getOrDefault(self % window, 'window', 10)
    print *, '<aqz22> [simpleFMClerkext] FM-window cycles set to:'
    print *, self % window

    ! NEW: Allocate fundamental eigenvector for scaling particles
    allocate(self % eigVec(self % N))   
    self % eigVec = ONE

    ! NEW: Allocate space and initialize the normalized/unnormalized matrix and startWgt
    allocate(self % tallyMatrix(self % window, self % N, self % N))
    allocate(self % matrix(self % N, self % N))
    allocate(self % startWgt(self % window, self % N))
    self % tallyMatrix = ZERO
    self % matrix = ZERO
    self % startWgt = ZERO

    ! Initialise response
    call self % resp % build(macroNuFission)

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual,'handleVirtual', .true.)

    ! NEW: extra prints:
    call dict % getOrDefault(self % doDebug, 'doDebug', .false.)
    if (self % doDebug) print *, '<aqz22> [simpleFMClerkext] DEBUG MESSAGES ENABLED' 
  
    ! NEW: force one:
    call dict % getOrDefault(self % forceOne, 'forceOne', .false.)
    if (self % forceOne) print *, '<aqz22> [simpleFMClerkext] FORCED HOMOGENEOUS EIGENVECTOR ENABLED' 

  

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(simpleFMClerkExt),intent(in)            :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleStart_CODE, closeCycle_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(simpleFMClerkExt), intent(in) :: self
    integer(shortInt)                :: S

    S = self % N * (self % N + 1)

  end function getSize

  !!
  !! Process start of the cycle
  !! Calculate starting weights in each bin and store them at memory location:
  !! self % getMemAddress() : self % getMemAddress() + N - 1
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleStart(self, start, mem)
    class(simpleFMClerkExt), intent(inout)       :: self
    class(particleDungeon), intent(in)           :: start
    type(scoreMemory), intent(inout)             :: mem
    integer(shortInt)                            :: idx, i

    if (self % doDebug) print *, '<aqz22> [simpleFMClerkext] cycle start'

    ! NEW: pop the first element of the queues and create new element for the next cycle (FiLo)
    do i = 2, self % window
      self % startWgt(i - 1, : ) = self % startWgt(i, : )
      self % tallyMatrix(i - 1, : , : ) = self % tallyMatrix(i, : , : )
    end do
    
    self % startWgt(self % window, : ) = ZERO
    self % tallyMatrix(self % window, : , : ) = ZERO


    ! NEW: Loop through a population, calculate starting weight in each bin, add to latest in moving window
    do i = 1, start % popSize()

      associate (state => start % get(i))

        idx = self % map % map(state)
        if (idx == 0) cycle
        call mem % score(state % wgt, self % getMemAddress() + idx - 1)
        self % startWgt(self % window, idx) = self % startWgt(self % window, idx) + state % wgt
      end associate

    end do

    ! NEW: print if necessary
    if (self % doDebug) then
      print *, '<aqz22> [simpleFMClerkext] origin bin weights:'
      do i = 1, self % window
        print *, self % startWgt(i, :)
      end do
    end if



  end subroutine reportCycleStart

  !!
  !! Process incoming collision report
  !!
  !! Calculate matrix elements and store them at memory location:
  !! self % getMemAddress() + N : self % getMemAddress() + N*(1 + N)
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(simpleFMClerkExt), intent(inout)  :: self
    class(particle), intent(in)          :: p
    class(nuclearDatabase),intent(inout) :: xsData
    type(scoreMemory), intent(inout)     :: mem
    logical(defBool), intent(in)         :: virtual
    class(neutronMaterial), pointer      :: mat
    type(particleState)                  :: state
    integer(shortInt)                    :: sIdx, cIdx
    integer(longInt)                     :: addr
    real(defReal)                        :: score, flux
    character(100), parameter :: Here = 'reportInColl simpleFMClerkExt_class.f90'

    ! Return if collision is virtual but virtual collision handling is off
    if ((.not. self % handleVirtual) .and. virtual) return

    ! Ensure we're not in void (could happen when scoring virtual collisions)
    if (p % matIdx() == VOID_MAT) return

    ! Get material pointer
    mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))
    if (.not.associated(mat)) then
      call fatalError(Here,'Unrecognised type of material was retrived from nuclearDatabase')
    end if

    ! Return if material is not fissile
    if (.not. mat % isFissile()) return

    ! Calculate flux with the right cross section according to virtual collision handling
    if (self % handleVirtual) then
      flux = p % w / xsData % getTrackingXS(p, p % matIdx(), TRACKING_XS)
    else
      flux = p % w / xsData % getTotalMatXS(p, p % matIdx())
    end if

    ! Find starting index in the map
    sIdx = self % map % map(p % preHistory)

    ! Find collision index in the map
    state = p
    cIdx = self % map % map(state)

    ! Defend against invalid collision or starting bin
    if (cIdx == 0 .or. sIdx == 0) return

    ! Calculate fission neutron production
    score = self % resp % get(p, xsData) * flux

    ! Score element of the matrix
    ! Note that the matrix memory location starts from memAddress + N
    addr = self % getMemAddress() + sIdx * self % N + cIdx - 1
    call mem % score(score, addr)

    ! NEW: Score to non-normalized matrix
    self % tallyMatrix(self % window, cIdx, sIdx) = self % tallyMatrix(self % window, cIdx, sIdx) + score
    

  end subroutine reportInColl

  !!
  !! Process cycle end
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine closeCycle(self, end, mem)
    class(simpleFMClerkExt), intent(inout) :: self
    class(particleDungeon), intent(in)  :: end
    type(scoreMemory), intent(inout)    :: mem
    integer(shortInt)                   :: i, j
    integer(longInt)                    :: addrFM
    real(defReal)                       :: normFactor


    if (mem % lastCycle()) then

      if (self % doDebug) print *, '<aqz22> [simpleFMClerkext] cycle end'

      ! Set address to the start of Fission Matrix
      ! Decrease by 1 to get correct address on the first iteration of the loop
      addrFM  = self % getMemAddress() + self % N - 1

      ! Normalise and accumulate estimates
      do i = 1, self % N

        ! Calculate normalisation factor
        normFactor = mem % getScore(self % getMemAddress() + i - 1)
        if (normFactor /= ZERO) normFactor = ONE / normFactor

        do j = 1, self % N
          ! Normalise FM column
          addrFM = addrFM + 1
          call mem % closeBin(normFactor, addrFM)
        end do

      end do






      ! New: normalize every cycle
      if (self % doDebug) print *, '<aqz22> [simpleFMClerkext] normalise FM factors:'

      do j = 1, self % N
        ! Calculate normalisation factor
        normFactor = sum(self % startWgt( : , j))
        if (normFactor /= ZERO) normFactor = ONE / normFactor
        
        if (self % doDebug) write(*,'(F18.15)', advance='no') normFactor
        
        do i = 1, self % N
          self % matrix(i,j) = sum(self % tallyMatrix( : , i,j)) * normFactor
        end do
      end do

      if (self % doDebug) print *, ''

      ! NEW: Obtain the fission matrix eigenvector
      call self % solve(mem)
    end if

  end subroutine closeCycle


  !!
  !! NEW:
  !! Solve the fission matrix eigenvalue problem
  !! by power iteration
  !!
  subroutine solve(self, mem)
    class(simpleFMClerkExt), intent(inout)       :: self
    real(defReal), dimension(:), allocatable     :: b
    real(defReal)                                :: tol, err
    integer(shortInt)                            :: it, i, j, itMax
    type(scoreMemory), intent(inout)             :: mem ! NEW: memory
    real(defReal), dimension(:,:), allocatable   :: totMat ! NEW
    
    if (self % doDebug) then
      
      allocate(totMat(self % N, self % N))
      do j = 1, self % N
        do i = 1, self % N
          totMat(i, j) = sum(self % tallyMatrix(:, i, j))
        end do
      end do

      print *, '<aqz22> [simpleFMClerkext] tally matrix:'
      print *, totMat
      print *, '<aqz22> [simpleFMClerkext] FM matrix:'
      print *, self % matrix

      print *, '<aqz22> [simpleFMClerkext] Start EV-solve...'
    end if


    tol = 1.0E-7
    err = ONE
    it = 0
    itMax = 10000
    allocate(b(self % N))
    self % eigVec = ONE


    do it = 1, itMax 

      b = self % eigVec 
      self % eigVec = ZERO

      ! Matrix-vector multiply
      do i = 1, self % N
        do j = 1, self % N
          self % eigVec(i) = self % eigVec(i) + self % matrix(i,j) * b(j)
        end do
      end do

      ! Normalise appropriately
      self % eigVec = self % eigVec / norm2(self % eigVec)
      
      err = norm2(self % eigVec - b) / norm2(b)
      if (err < tol .and. it > 200) exit


    end do

    if (it >= itMax) print *,'WARNING: FM iterations did not finish <aqz22> [simpleFMClerkext]'
  
    self % eigVec = self % eigVec / sum(self % eigVec)

    if (self % doDebug) then
      print *,'<aqz22> [simpleFMClerkext] EV finished. Iterations: '
      print *, it
    end if

    if (self % forceOne) self % eigVec = ONE

    print *,'FM Eigenvector:'
    print *, self % eigVec


  end subroutine solve



  !!
  !! Return result from the clerk for interaction with Physics Package
  !!  Returns FMresult defined in this module
  !!   If res is already allocated to a FM of fitting size it reuses already allocated space
  !!    This should improve performance when updating estimate of FM each cycle
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(simpleFMClerkExt), intent(in)               :: self
    class(tallyResult),allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                  :: mem
    integer(shortInt)                              :: i, j
    integer(longInt)                               :: addr
    real(defReal)                                  :: val, STD

    !! Allocate result to FMmatrix
    !! Do not deallocate if already allocated to FMmatrix
    !! Its not to nice -> clean up
    !if (allocated(res)) then
!
    !  select type(res)
    !    class is (FMmatrix)
    !      ! Do nothing
    !    class default
    !      ! Reallocate
    !      deallocate(res)
    !      allocate( FMmatrix :: res)
    !  end select
!
    !else
    !  allocate( FMmatrix :: res)
!
    !end if
!
    !! Load data into the FM
    !select type(res)
    !  class is(FMmatrix)
    !    ! Check size and reallocate space if needed
    !    ! This is horrible. Hove no time to polish. Blame me (MAK)
    !    if (allocated(res % FM)) then
!
    !      if (any(shape(res % FM) /= [self % N, self % N, 2])) then
    !        deallocate(res % FM)
    !        allocate(res % FM(self % N, self % N, 2))
    !      end if
!
    !    else
    !      allocate(res % FM(self % N, self % N, 2))
    !    end if
!
    !    ! Set size of the FM
    !    res % N = self % N
!
    !    ! Load entries
    !    addr = self % getMemAddress() + self % N - 1
    !    do i = 1, self % N
    !      do j = 1, self % N
    !        addr = addr + 1
    !        call mem % getResult(val, STD, addr)
    !        res % FM(j, i, 1) = val
    !        res % FM(j, i, 2) = STD
    !      end do
    !    end do
!
    !end select
!


    !! NEW: also return the resulting FM-ev
        ! Allocate result to FMeigen
    ! Do not deallocate if already allocated to FMeigen
    ! Its not to nice -> clean up
    if (allocated(res)) then
      select type(res)
        class is (FMresult)
          ! Do nothing

        class default ! Reallocate
          deallocate(res)
          allocate( FMresult :: res)
     end select

    else
      allocate( FMresult :: res)

    end if

    ! Load data inti the FM
    select type(res)
      class is(FMresult)
        ! Check size and reallocate space if needed
        ! This is horrible. Hove no time to polish. Blame me (MAK)
        if (allocated(res % eigVec)) then
          if (any(shape(res % eigVec) /= self % N)) then
            deallocate(res % eigVec)
            allocate(res % eigVec(self % N))
          end if
        else
          allocate(res % eigVec(self % N))
        end if

        ! Set size of the FM
        res % N = self % N

        ! Load entries
        do i = 1,self % N
          res % eigVec(i) = self % eigVec(i)
        end do

    end select

  end subroutine getResult



  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(simpleFMClerkExt), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem

    call statusMsg('simpleFMClerkExt does not support display yet')

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(simpleFMClerkExt), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    integer(shortInt)                :: i
    integer(longInt)                 :: addr
    real(defReal)                    :: val, std
    character(nameLen)               :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print map information
    call self % map % print(outFile)

    ! Print fission matrix
    name = 'FM'
    addr = self % getMemAddress() + self % N - 1

    call outFile % startArray(name, [self % N, self % N])

    do i = 1, self % N * self % N
      addr = addr + 1
      call mem % getResult(val, std, addr)
      call outFile % addResult(val, std)
    end do
    call outFile % endArray()


    ! NEW: print eigenvector
    name = 'eig'
    call outFile % startArray(name, [self % N, 1])
    do i = 1,self % N
      val = self % eigVec(i)
      call outFile % addResult(val, ZERO)
    end do
    call outFile % endArray()


    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  !! See tallyClerk_inter for details
  !!
  elemental subroutine kill(self)
    class(simpleFMClerkExt), intent(inout) :: self

    ! Call superclass
    call kill_super(self)

    if (allocated(self % map)) deallocate(self % map)

    ! NEW: deallocate matrices
    if (allocated(self % tallyMatrix)) deallocate(self % tallyMatrix)
    if (allocated(self % matrix)) deallocate(self % matrix)
    if (allocated(self % startWgt)) deallocate(self % startWgt)

    self % N = 0
    self % handleVirtual = .true.

    call self % resp % kill()

  end subroutine kill

end module simpleFMClerkExt_class
