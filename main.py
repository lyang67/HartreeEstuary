# This is a sample Python script.
import math

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
gravity = 9.8

class GridPoint:
    def __init__(self, d, c, v):
        self.depth = d
        self.celerity = c
        self.velocity = v


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def set_init_conditions():
    leftBoundaryGridPoint = GridPoint(1.4, math.sqrt(gravity * 1.4), 0.1)
    site2CGridpoint = GridPoint(1.3, math.sqrt(gravity * 1.3), 0.1)
    site22GridPoint = GridPoint(1.2, math.sqrt(gravity * 1.2), 0.2)
    site21GridPoint = GridPoint(1.1, math.sqrt(gravity * 1.1), 0.3)
    site2BGridPoint = GridPoint(1.0, math.sqrt(gravity * 1.0), 0.4)
    site2GridPoint = GridPoint(0.9, math.sqrt(gravity * 0.9), 0.5)
    rightBoundaryGridPoint = GridPoint(0.8, math.sqrt(gravity * 0.8), 0.1)

    initGridPoints = [leftBoundaryGridPoint, site2CGridpoint, site22GridPoint, site21GridPoint, site2BGridPoint,
                      site2GridPoint, rightBoundaryGridPoint]
    print(initGridPoints)

    allGridPoints = [initGridPoints]
    return allGridPoints

# only get the left characteristic
def calculateRightBoundaryConditions(time, gridPointA, gridpointB, dtDx):
    if (time <= 0):
        return

    #TODO make a proper function for tidal depths according to time, just use a rough placeholder function now
    rightBoundaryDepth = math.sin(4.7 + 0.25 * time) + 1.2
    rightBoundaryCelerity = math.sqrt(gravity * rightBoundaryDepth)

    leftGridPoint = calculateLeft(dtDx, gridPointA, gridpointB)
    velocityLeft = leftGridPoint.velocity

    celerityLeft = leftGridPoint.celerity

    rightBoundaryVelocity = velocityLeft + 2*(celerityLeft - rightBoundaryCelerity)

    rightBoundaryGridpoint = GridPoint(rightBoundaryDepth, rightBoundaryCelerity, rightBoundaryVelocity)
    return rightBoundaryGridpoint

def calculateLeftBoundaryConditions(time, gridPointB, gridpointC, dtDx):
    #TODO: assuming we know depth at this left boundary but unsure what we know here, check with H Chanson
    if (time <= 0):
        return

    #TODO make a proper function for tidal depths according to time, just use a rough placeholder function now
    leftBoundaryDepth = math.sin(0.2 * time) + 1.4
    leftBoundaryCelerity = math.sqrt(gravity * leftBoundaryDepth)

    rightGridPoint = calculateRight(dtDx, gridPointB, gridpointC)
    velocityright = rightGridPoint.velocity

    celerityright = rightGridPoint.celerity

    leftBoundaryVelocity = velocityright + 2*(celerityright - leftBoundaryCelerity)

    leftBoundaryGridpoint = GridPoint(leftBoundaryDepth, leftBoundaryCelerity, leftBoundaryVelocity)
    return leftBoundaryGridpoint

def calculateLeft(dtDx, gridPointA, gridpointB):
    velocityLeft = (gridpointB.velocity +
                    dtDx*(gridpointB.celerity * gridPointA.velocity - gridPointA.celerity * gridpointB.velocity)) \
                   / (1 + dtDx*(gridpointB.velocity +gridpointB.celerity - gridPointA.velocity - gridPointA.celerity))

    celerityLeft = ((gridpointB.celerity) + velocityLeft*(gridPointA.celerity - gridpointB.celerity)) \
                   / (1 + dtDx * (gridpointB.celerity - gridPointA.celerity))

    depthLeft = (celerityLeft ** 2) / gravity

    gridPointLeft = GridPoint(depthLeft, celerityLeft, velocityLeft)

    return gridPointLeft


def calculateRight(dtDx, gridpointB, gridpointC):
    velocityRight = (gridpointB.velocity +
                    dtDx * (gridpointB.celerity * gridpointC.velocity - gridpointC.celerity * gridpointB.velocity)) \
                   / (1 + dtDx * (
                gridpointC.velocity + gridpointC.celerity - gridpointB.velocity - gridpointB.celerity))

    celerityRight = ((gridpointB.celerity) + velocityRight * (gridpointB.celerity - gridpointC.celerity)) \
                   / (1 + dtDx * (gridpointB.celerity - gridpointC.celerity))

    depthRight = (celerityRight ** 2) / gravity

    gridPointRight = GridPoint(depthRight, celerityRight, velocityRight)

    return gridPointRight

def populateGrid(mainGrid, timestepSize, xDistanceSize):
    # calculate the number of timesteps we have in a 24 hour period
    numTimesteps = (24 * 3600)/timestepSize

    # we have 5 known stations, and 2 boundary conditions, so 5 + 2 grids points in the horizontal direction
    numHorizontalGridPoints = 7

    # 0 based indexing for the grid
    currentTimestep = 1
    currentGridPoint = 1

    while (currentTimestep <= numTimesteps):
        # append a new empty list of gridpoints to the current time step
        mainGrid.append([])

        # start populating the current list of gridpoints at this timestep
        leftBoundaryConditions = calculateLeftBoundaryConditions\
            (timestepSize * currentTimestep,
             mainGrid[currentTimestep - 1][currentGridPoint -1],
             mainGrid[currentTimestep - 1][currentGridPoint],
             timestepSize/xDistanceSize)

        mainGrid[currentTimestep].append(leftBoundaryConditions)


        #calculate conditions at each grid point between boundary
        #gridpoints 1 and 7 are the boundary conditions points
        while (currentGridPoint < numHorizontalGridPoints - 1):
            #calculate flow conditions at this grid point
            #get grid points at (x-1), (t-1) and x, t-1 and x+1, t-1, name them A, B C as in textbook convention
            #shouldn't have to worry about indexing beyond bounds of list as we are starting at 1+ our boundary conditions
            gridPointA = mainGrid[currentTimestep - 1][currentGridPoint - 1]
            gridPointC = mainGrid[currentTimestep - 1][currentGridPoint + 1]
            gridPointB = mainGrid[currentTimestep - 1][currentGridPoint]

            left = calculateLeft(timestepSize/xDistanceSize, gridPointA, gridPointB)
            right = calculateRight(timestepSize/xDistanceSize, gridPointB, gridPointC)

            # ignore friction factor for now
            currentVelocity = (left.velocity + right.velocity + 2*(left.celerity - right.celerity)) / 2
            currentCelerity = (left.velocity - right.velocity + 2*(left.celerity + right.celerity)) / 4
            currentDepth = (currentCelerity ** 2) / gravity

            mGridPoint = GridPoint(currentDepth, currentCelerity, currentVelocity)
            mainGrid[currentTimestep].append(mGridPoint)
            currentGridPoint += 1

        # we have incremented the index to the boundary condition index before exiting loop, so point "A" is at gridpoint -1 \
        # and point "B" is at current grid point
        rightBoundaryConditions = calculateRightBoundaryConditions(
            timestepSize * currentTimestep,
            mainGrid[currentTimestep - 1][currentGridPoint - 1],
            mainGrid[currentTimestep - 1][currentGridPoint],
            timestepSize / xDistanceSize
        )
        mainGrid[currentTimestep].append(rightBoundaryConditions)
        currentGridPoint = 1
        currentTimestep += 1

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #TODO: the sites don't look 250m apart, but figure this out later
    gridSizeX = 250
    #every 1 hour (in units seconds) for now, so we don't have to wait forever for the computation to finish while testing
    gridSizeY = 3600
    hartTreeGrid = set_init_conditions()
    populateGrid(hartTreeGrid,gridSizeY, gridSizeX)

    timestepCounter = 0
    for timestep in hartTreeGrid:
        print('Timestep number ' + str(timestepCounter))
        gridCounter = 0
        for gridpoint in hartTreeGrid[timestepCounter]:
            print('Grid number ' + str(gridCounter))
            print('Celerity is ' + str(gridpoint.celerity))
            print('Velocity is ' + str(gridpoint.velocity))
            print('Depth is ' + str(gridpoint.depth))
            gridCounter += 1
        timestepCounter += 1

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
