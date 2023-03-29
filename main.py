# This is a sample Python script.
import math
import csv

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
gravity = 9.8

class GridPoint:
    def __init__(self, d, c, v, x=0, ca=0, pw=0):
        self.depth = d
        self.celerity = c
        self.velocity = v
        self.crossSectionArea = ca
        self.wettedPerimeter = pw
        self.xLocation = x

def set_init_conditions():
    initialDepth = 2.3
    initialCelerity = math.sqrt(gravity * initialDepth)

    initialCrossSection = initialDepth*50
    initialWettedPerimeter = 2*initialDepth + 50

    leftBoundaryGridPoint = GridPoint(initialDepth, initialCelerity, 0.1, 0, initialCrossSection, initialWettedPerimeter)
    site2CGridpoint = GridPoint(initialDepth, initialCelerity, 0.1, 300, initialCrossSection, initialWettedPerimeter)
    site22GridPoint = GridPoint(initialDepth, initialCelerity, 0.1, 600, initialCrossSection, initialWettedPerimeter)
    site21GridPoint = GridPoint(initialDepth, initialCelerity, 0.1, 900, initialCrossSection, initialWettedPerimeter)
    site2BGridPoint = GridPoint(initialDepth, initialCelerity, 0.1, 1200, initialCrossSection, initialWettedPerimeter)
    site2GridPoint = GridPoint(initialDepth, initialCelerity, 0.1, 1500, initialCrossSection, initialWettedPerimeter)
    rightBoundaryGridPoint = GridPoint(initialDepth, initialCelerity, 0.1, 1800, initialCrossSection, initialWettedPerimeter)

    initGridPoints = [leftBoundaryGridPoint, site2CGridpoint, site22GridPoint, site21GridPoint, site2BGridPoint,
                      site2GridPoint, rightBoundaryGridPoint]
    print(initGridPoints)

    allGridPoints = [initGridPoints]
    return allGridPoints

# only get the left characteristic
def calculateRightBoundaryConditions(time, gridPointA, gridpointB, dtDx):
    if (time <= 0):
        return

    #convert time to hours
    timeHours = time / 3600
    #TODO make a proper function for tidal depths according to time, just use a rough placeholder function now
    rightBoundaryDepth = math.cos(0.166 * math.pi * timeHours) + 1.3
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

    #left boundary velocity is 0 since there is no inflow
    leftBoundaryVelocity = 0

    # calculate right gridpoint
    rightGridPoint = calculateRight(dtDx, gridPointB, gridpointC)
    velocityright = rightGridPoint.velocity
    celerityright = rightGridPoint.celerity

    # calculate left celerity using backwards characteristic
    leftBoundaryCelerity = (leftBoundaryVelocity - velocityright + 2*celerityright) / 2
    leftBoundaryDepth = leftBoundaryCelerity**2 / gravity

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

def calculateWettedPerimeter(depth, xLocation):
    width = getBottomWidthAtLocation(xLocation)
    # rectangular channel only now
    return 2 * depth + width

def calculateHydraulicDiameter(depth, xLocation):
    # calculate stuff
    width = getBottomWidthAtLocation(xLocation)

    crossSectionArea = width * depth
    wettedPerimeter = calculateWettedPerimeter(depth, xLocation)

    return (4 * crossSectionArea) / wettedPerimeter

def calculateFrictionSlope(frictionFactor, interpolatedVelocity, hydraulicDiameter):
    frictionSlope = (frictionFactor * interpolatedVelocity**2)/(8 * gravity * (hydraulicDiameter / 4))

    return frictionSlope

def getBottomWidthAtLocation(xLocation):
    return 50

def populateGrid(mainGrid, timestepSize, xDistanceSize):
    # calculate the number of timesteps we have in a 24 hour period
    numTimesteps = (24 * 3600)/timestepSize
    dxDt = xDistanceSize / timestepSize

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

        if CheckCourantCondition(leftBoundaryConditions, dxDt) is not True:
            print('Courant condition failed at left boundary')
            return False

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
            xLocation = gridPointB.xLocation

            left = calculateLeft(timestepSize/xDistanceSize, gridPointA, gridPointB)
            right = calculateRight(timestepSize/xDistanceSize, gridPointB, gridPointC)

            # TODO update the left and right x location
            frictionFactor = 0.025
            hydraulicDiameterLeft = calculateHydraulicDiameter(left.depth, xLocation)
            frictionSlopeLeft = calculateFrictionSlope(frictionFactor, left.velocity, hydraulicDiameterLeft)

            hydraulicDiameterRight = calculateHydraulicDiameter(right.depth, xLocation)
            frictionSlopeRight = calculateFrictionSlope(frictionFactor, right.velocity, hydraulicDiameterRight)

            #alternative versions of velocity and celerity calculation with friction
            currentVelocity = (left.velocity + right.velocity + 2 * (left.celerity - right.celerity)) / 2 \
                              - gravity * timestepSize * (frictionSlopeLeft + frictionSlopeRight)
            currentCelerity = (left.velocity - right.velocity + 2 * (left.celerity + right.celerity)) / 4 \
                              - gravity * timestepSize * (frictionSlopeLeft - frictionSlopeRight)

            # ignore friction factor for now
            #currentVelocity = (left.velocity + right.velocity + 2*(left.celerity - right.celerity)) / 2
            #currentCelerity = (left.velocity - right.velocity + 2*(left.celerity + right.celerity)) / 4
            currentDepth = (currentCelerity ** 2) / gravity
            #same xLocation
            mGridPoint = GridPoint(currentDepth, currentCelerity, currentVelocity, xLocation)

            if CheckCourantCondition(mGridPoint, dxDt) is not True:
                print('Courant condition failed at in the middle')
                return False

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

        if CheckCourantCondition(mGridPoint, dxDt) is not True:
            print('Courant condition failed at right boundary')
            return False

        mainGrid[currentTimestep].append(rightBoundaryConditions)
        currentGridPoint = 1
        currentTimestep += 1
    return True

def CheckCourantCondition(gridpoint, dxDt):
    courantCoundition1 = abs(gridpoint.velocity + gridpoint.celerity) / dxDt
    courantCoundition2 = abs(gridpoint.velocity - gridpoint.celerity) / dxDt

    print('courant condition check left is ' + str(courantCoundition1) + ' and right is ' + str(courantCoundition2))
    if courantCoundition1 > 1 or courantCoundition2 > 1:
        return False

    return True

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #TODO: the sites don't look 250m apart, but figure this out later
    gridSizeX = 300
    #set our timestep to 30s
    gridSizeY = 30
    hartTreeGrid = set_init_conditions()
    populateGridResult = populateGrid(hartTreeGrid,gridSizeY, gridSizeX)

    if populateGridResult is not True:
        print('Courant condition check failed, check output')

    with open('results.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        headingGrids = [' ', 'Left Boundary', ' ',
                        '', 'Gridpoint 1', '',
                        '', 'Gridpoint 2', '',
                        '', 'Gridpoint 3', '',
                        '', 'Gridpoint 4', '',
                        '', 'Gridpoint 5', '',
                        '', 'Right Boundary', '']
        writer.writerow(headingGrids)

        headingTypes = ['Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)',
                        'Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)',
                        'Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)',
                        'Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)',
                        'Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)',
                        'Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)',
                        'Depth (m)', 'Celerity(m/s)', 'Velocity (m/s)']
        writer.writerow(headingTypes)

        timestepCounter = 0

        for timestep in hartTreeGrid:
            print('Timestep number ' + str(timestepCounter))
            gridCounter = 0
            currentTimestepValues = []
            for gridpoint in hartTreeGrid[timestepCounter]:
                currentTimestepValues.append(gridpoint.depth)
                currentTimestepValues.append(gridpoint.celerity)
                currentTimestepValues.append(gridpoint.velocity)
                #print('Grid number ' + str(gridCounter) + ' Celerity is ' + str(gridpoint.celerity)
                      #+ ' Depth is ' + str(gridpoint.depth) + ' Velocity is ' + str(gridpoint.velocity))
                gridCounter += 1
            writer.writerow(currentTimestepValues)
            timestepCounter += 1

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
