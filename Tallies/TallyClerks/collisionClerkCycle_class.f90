module collisionClerkCycle_class

  use numPrecision
  use tallyCodes
  use universalVariables
  use genericProcedures,          only : fatalError
  use display_func,               only : statusMsg
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill
  
  ! NEW: use particle dungeon
  use particleDungeon_class,      only : particleDungeon

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Responses
  use tallyResponseSlot_class,    only : tallyResponseSlot

  implicit none
  private

  !!
  !! Collision estimator of reaction rates
  !! Calculates flux weighted integral from collisions
  !!
  !! Private Members:
  !!   filter   -> Space to store tally Filter
  !!   map      -> Space to store tally Map
  !!   response -> Array of responses
  !!   width    -> Number of responses (# of result bins for each map position)
  !!   handleVirtual -> score on virtual collisions (due to TMS or delta tracking)
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myCollisionClerk {
  !!   type collisionClerk;
  !!   # handleVirtual 0; # default is 1   
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #)
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: collisionClerkCycle
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response

    ! NEW: track cycles
    integer(shortInt)             :: maxCycles = 1    !! Number of tally cycles
    integer(shortInt)             :: currentCycle = 0 !! track current cycle

    ! NEW: the memory size of one cycle
    integer(longInt)   :: cycleSize = 0

    ! Useful data
    integer(shortInt)  :: width = 0

    ! Settings
    logical(defBool)   :: handleVirtual = .true.

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! NEW: File reports and check status -> run-time procedures
    procedure  :: reportCycleEnd

    ! File reports and check status -> run-time procedures
    procedure  :: reportInColl

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type collisionClerkCycle

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(collisionClerkCycle), intent(inout)        :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen),dimension(:),allocatable :: responseNames
    integer(shortInt)                           :: i

    ! Assign name
    call self % setName(name)

    ! Load filetr
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

    ! NEW: Read number of cycles for which to track entropy
    call dict % get(self % maxCycles, 'cycles')

    ! Get names of response dictionaries
    call dict % get(responseNames,'response')

    ! Load responses
    allocate(self % response(size(responseNames)))
    do i=1, size(responseNames)
      call self % response(i) % init(dict % getDictPtr( responseNames(i) ))
    end do

    ! Set width
    self % width = size(responseNames)

    ! NEW: store size of one cycle
    self % cycleSize = size(self % response)
    if(allocated(self % map)) self % cycleSize = self % cycleSize * self % map % bins(0) 

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual,'handleVirtual', .true.)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(collisioNClerkCycle), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if (allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    ! Kill and deallocate responses
    if (allocated(self % response)) then
      deallocate(self % response)
    end if

    self % width   = 0
    self % handleVirtual = .true.

    ! NEW: variables
    self % currentCycle = 0
    self % maxCycles = 1
    self % cycleSize = 0

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(collisionClerkCycle),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    ! NEW: can also return the cycle
    validCodes = [inColl_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(collisionClerkCycle), intent(in) :: self
    integer(shortInt)                      :: S

    S = size(self % response)
    if(allocated(self % map)) S = S * self % map % bins(0) 

    ! NEW: size for all cycles
    S = S * self % maxCycles

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(collisionClerkCycle), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    logical(defBool), intent(in)          :: virtual
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx, i
    integer(longInt)                      :: addr
    real(defReal)                         :: scoreVal, flux
    character(100), parameter :: Here = 'reportInColl (collisionClerkCycle_class.f90)'

    ! Return if collision is virtual but virtual collision handling is off
    if ((.not. self % handleVirtual) .and. virtual) return

    ! NEW: return if maximum cycle is reached
    if (self % currentCycle >= self % maxCycles) return

    ! Get current particle state
    state = p

    ! Check if within filter
    if (allocated(self % filter)) then
      if (self % filter % isFail(state)) return
    end if

    ! Find bin index
    if (allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate flux with the right cross section according to virtual collision handling
    if (self % handleVirtual) then
      flux = p % w / xsData % getTrackingXS(p, p % matIdx(), TRACKING_XS)
    else
      flux = p % w / xsData % getTotalMatXS(p, p % matIdx())
    end if

    ! Calculate bin address
    addr = self % getMemAddress() + self % width * (binIdx - 1)  - 1

    ! NEW: increase address to match the current cycle
    addr = addr + self % currentCycle * self % cycleSize

    ! Append all bins
    do i = 1, self % width
      scoreVal = self % response(i) % get(p, xsData) * flux
      call mem % score(scoreVal, addr + i)

    end do

  end subroutine reportInColl


  !!
  !! NEW: reportCycleEnd, increments the cycle number
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(collisionClerkCycle), intent(inout) :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem

    self % currentCycle = self % currentCycle + 1

  end subroutine reportCycleEnd


  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(collisionClerkCycle), intent(in)  :: self
    type(scoreMemory), intent(in)      :: mem

    call statusMsg('collisionClerkCycle does not support display yet')

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(collisionClerkCycle), intent(in)          :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information
    if (allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array
    ! NEW: added , self % maxCycles
    if (allocated(self % map)) then
      resArrayShape = [size(self % response), self % map % binArrayShape(), self % maxCycles]
    else
      resArrayShape = [size(self % response), self % maxCycles]
    end if

    ! Start array
    name ='Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i = 1, product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val,std)

    end do

    call outFile % endArray()
    call outFile % endBlock()

  end subroutine print

end module collisionClerkCycle_class
