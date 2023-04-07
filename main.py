# This is a sample Python script.
import math
import csv

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
gravity = 9.8
frictionFactor = 0.02

# we have 5 known stations, and 2 boundary conditions, so 5 + 2 grids points in the horizontal direction
#numHorizontalGridPoints = 10
numHorizontalGridPoints = 121
#original 300 dx, 30 dt
horizontalGridSize = 0.1
verticalGridSize = 0.05
#channelWidth = 12
channelWidth = 0.5
#initialDepth = 6
#initialDepth = 2
initialDepth = 0.175778027747
initialQ = 0.03448
#hoursToRun = 1
hoursToRun = 14/60
depthData = []
withFriction = True
diagnosticPrint = False
compareAnalytical = False
#big ugly switch to turn off everything in main
testOther = False

class GridPoint:
    def __init__(self, d, c, v, sf, x):
        self.depth = d
        self.celerity = c
        self.velocity = v
        self.frictionSlope = sf
        self.xLocation = x

def populateDepthBoundary():
    with open('DepthData.csv', mode='r') as file:
        # reading the CSV file
        depthFile = csv.reader(file)

        # displaying the contents of the CSV file
        for lines in depthFile:
            depthData.append(float(lines[0]))
        print (depthData)

def set_init_conditions():
    initialCelerity = math.sqrt(gravity * initialDepth)
    #initialVelocity = 0
    initialCrossSection = initialDepth*channelWidth
    initialVelocity = -initialQ /initialCrossSection
    initialWettedPerimeter = 2*initialDepth + channelWidth
    initialHydraulicDiameter = (4*initialCrossSection) / initialWettedPerimeter
    initialFrictionSlope = calculateFrictionSlope(frictionFactor, initialVelocity, initialHydraulicDiameter)

    gridCounter = 0
    initialGridPoints = []
    while gridCounter < numHorizontalGridPoints:
        initialGridPoint = GridPoint(initialDepth, initialCelerity, initialVelocity, initialFrictionSlope,
                                     gridCounter * horizontalGridSize)
        initialGridPoints.append(initialGridPoint)
        gridCounter += 1

    allGridPoints = [initialGridPoints]
    return allGridPoints

# only get the left characteristic
def calculateLeftBoundaryConditions(time, gridPointB, gridpointC, dt, dx):
    if (time <= 0):
        return

    dtDx = dt / dx

    #convert time to hours
    timeHours = time / 3600
    #TODO make a proper function for tidal depths according to time, just use a rough placeholder function now
    #rightBoundaryDepth = math.cos(0.166 * math.pi * timeHours) + 6
    #leftBoundaryDepth = 6-(time/3600)*0.5
    depthIndex = int(time/verticalGridSize)
    leftBoundaryDepth = depthData[depthIndex]
    leftBoundaryCelerity = math.sqrt(gravity * leftBoundaryDepth)

    rightGridPoint = calculateRight(dtDx, gridPointB, gridpointC)
    velocityRight = rightGridPoint.velocity

    celerityRight = rightGridPoint.celerity

    if withFriction:
        rightBoundaryVelocity = velocityRight - 2*celerityRight + 2*leftBoundaryCelerity - \
                                gravity*dt * gridPointB.frictionSlope
    else:
        rightBoundaryVelocity = velocityRight - 2*celerityRight + 2*leftBoundaryCelerity

    currentHydraulicDiameter = calculateHydraulicDiameter(leftBoundaryDepth, gridPointB.xLocation)
    currentFrictionSlope = calculateFrictionSlope(frictionFactor, rightBoundaryVelocity, currentHydraulicDiameter)
    leftBoundaryGridpoint = GridPoint(leftBoundaryDepth, leftBoundaryCelerity, rightBoundaryVelocity,
                                       currentFrictionSlope, gridPointB.xLocation)
    return leftBoundaryGridpoint

def calculateRightBoundaryConditions(time, gridPointA, gridpointB, dt, dx):
    if (time <= 0):
        return
    dtDx = dt / dx
    #right boundary velocity is 0 since there is no inflow
    #rightBoundaryVelocity = 0

    if time < 120:
        rightBoundaryVelocity = (1/120) * time
    else:
        rightBoundaryVelocity = 1

    #boundary condiiton
    rightBoundaryDepth = initialDepth
    rightBoundaryCelerity = math.sqrt(gravity * rightBoundaryDepth)

    # calculate right gridpoint
    leftGridPoint = calculateLeft(dtDx, gridPointA, gridpointB)
    velocityLeft= leftGridPoint.velocity
    celerityLeft = leftGridPoint.celerity

    if withFriction:
        # calculate left celerity using backwards characteristic
        #rightBoundaryCelerity = (velocityLeft + 2*celerityLeft -
                                #gravity * dt * (1 * gridPointA.frictionSlope)) / 2

        rightBoundaryVelocity = velocityLeft +2*celerityLeft - 2*gridpointB.celerity - gravity * dt * (gridpointB.frictionSlope)
    else:
        rightBoundaryCelerity = (velocityLeft + 2*celerityLeft) / 2

    #rightBoundaryDepth = rightBoundaryCelerity**2 / gravity

    currentHydraulicDiameter = calculateHydraulicDiameter(rightBoundaryDepth, gridPointA.xLocation)
    currentFrictionSlope = calculateFrictionSlope(frictionFactor, rightBoundaryVelocity, currentHydraulicDiameter)
    rightBoundaryGridpoint = GridPoint(rightBoundaryDepth, rightBoundaryCelerity, rightBoundaryVelocity,
                                      currentFrictionSlope, gridPointA.xLocation)
    return rightBoundaryGridpoint

def calculateLeft(dtDx, gridPointA, gridpointB):
    velocityLeft = (gridpointB.velocity +
                    dtDx*(gridpointB.celerity * gridPointA.velocity - gridPointA.celerity * gridpointB.velocity)) \
                   / (1 + dtDx*(gridpointB.velocity +gridpointB.celerity - gridPointA.velocity - gridPointA.celerity))

    celerityLeft = ((gridpointB.celerity) + dtDx* velocityLeft*(gridPointA.celerity - gridpointB.celerity)) \
                   / (1 + dtDx * (gridpointB.celerity - gridPointA.celerity))

    depthLeft = (celerityLeft ** 2) / gravity

    #todo for completeness calculate xL

    #don't care about friction slope or location on the characteristic for now
    gridPointLeft = GridPoint(depthLeft, celerityLeft, velocityLeft, 0, 0)

    return gridPointLeft


def calculateRight(dtDx, gridpointB, gridpointC):
    velocityRight = (gridpointB.velocity +
                    dtDx * (gridpointB.celerity * gridpointC.velocity - gridpointC.celerity * gridpointB.velocity)) \
                   / (1 + dtDx * (
                gridpointC.velocity + gridpointB.celerity - gridpointB.velocity - gridpointC.celerity))

    celerityRight = ((gridpointB.celerity) + dtDx * velocityRight * (gridpointB.celerity - gridpointC.celerity)) \
                   / (1 + dtDx * (gridpointB.celerity - gridpointC.celerity))

    depthRight = (celerityRight ** 2) / gravity

    gridPointRight = GridPoint(depthRight, celerityRight, velocityRight, 0, 0)

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
    #frictionSlope = (frictionFactor * interpolatedVelocity * abs(interpolatedVelocity))/(8 * gravity * (hydraulicDiameter / 4))
    frictionSlope = (frictionFactor * interpolatedVelocity**2) / (
                8 * gravity * (hydraulicDiameter / 4))

    return frictionSlope

#TODO replace with varied bottom width function
def getBottomWidthAtLocation(xLocation):
    return channelWidth

def populateGrid(mainGrid, timestepSize, xDistanceSize):
    # calculate the number of timesteps we have in a 24 hour period
    numTimesteps = (hoursToRun * 3600)/timestepSize
    dxDt = xDistanceSize / timestepSize


    # 0 based indexing for the grid
    currentTimestep = 1
    currentGridPoint = 1

    while (currentTimestep <= numTimesteps):

        print('Timestep number ' + str(currentTimestep))
        # append a new empty list of gridpoints to the current time step
        mainGrid.append([])

        # start populating the current list of gridpoints at this timestep
        leftBoundaryConditions = calculateLeftBoundaryConditions\
            (timestepSize * currentTimestep,
             mainGrid[currentTimestep - 1][currentGridPoint -1],
             mainGrid[currentTimestep - 1][currentGridPoint],
             timestepSize, xDistanceSize)

        if CheckCourantCondition(leftBoundaryConditions, dxDt) is not True:
            print('Courant condition failed at left boundary')
            return False

        mainGrid[currentTimestep].append(leftBoundaryConditions)


        #calculate conditions at each grid point between boundary
        #gridpoints 1 and 7 are the boundary conditions points
        while (currentGridPoint < numHorizontalGridPoints - 1):
            print('current grid point is ' + str(currentGridPoint))
            #calculate flow conditions at this grid point
            #get grid points at (x-1), (t-1) and x, t-1 and x+1, t-1, name them A, B C as in textbook convention
            #shouldn't have to worry about indexing beyond bounds of list as we are starting at 1+ our boundary conditions
            gridPointA = mainGrid[currentTimestep - 1][currentGridPoint - 1]
            gridPointC = mainGrid[currentTimestep - 1][currentGridPoint + 1]
            gridPointB = mainGrid[currentTimestep - 1][currentGridPoint]
            xLocation = gridPointB.xLocation

            left = calculateLeft(timestepSize/xDistanceSize, gridPointA, gridPointB)
            right = calculateRight(timestepSize/xDistanceSize, gridPointB, gridPointC)

            if withFriction:
                # TODO update the left and right x location
                #alternative versions of velocity and celerity calculation with friction
                currentVelocity = (left.velocity + right.velocity + 2 * (left.celerity - right.celerity) \
                                   - gravity * timestepSize * gridPointB.frictionSlope) / 2
                currentCelerity = (left.velocity - right.velocity + 2 * (left.celerity + right.celerity))  / 4

            else:
                # ignore friction factor for now
                currentVelocity = (left.velocity + right.velocity + 2*(left.celerity - right.celerity)) / 2
                currentCelerity = (left.velocity - right.velocity + 2*(left.celerity + right.celerity)) / 4

            currentDepth = (currentCelerity ** 2) / gravity

            #same xLocation
            currentHydraulicDiameter = calculateHydraulicDiameter(currentDepth, gridPointB.xLocation)
            currentFrictionSlope = calculateFrictionSlope(frictionFactor, currentVelocity, currentHydraulicDiameter)
            mGridPoint = GridPoint(currentDepth, currentCelerity, currentVelocity, currentFrictionSlope,
                                   gridPointB.xLocation)

            print('velocity is ' + str(currentVelocity) + ' and celerity is ' + str(currentCelerity) +
                  ' and depth is ' +str(currentDepth) + ' and friction slope is ' + str(currentFrictionSlope))
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
            timestepSize, xDistanceSize
        )

        print('Right boundary velocity is ' + str(rightBoundaryConditions.velocity) + ' and celerity is ' + str(rightBoundaryConditions.celerity) +
              ' and depth is ' + str(rightBoundaryConditions.depth) + ' and friction slope is ' + str(rightBoundaryConditions.frictionSlope))

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

    if not testOther:

        populateDepthBoundary()
        #TODO: the sites don't look 250m apart, but figure this out later
        gridSizeX = horizontalGridSize
        #set our timestep to 30s
        gridSizeY = verticalGridSize
        hartTreeGrid = set_init_conditions()
        populateGridResult = populateGrid(hartTreeGrid,gridSizeY, gridSizeX)

        if populateGridResult is not True:
            print('Courant condition check failed, check output')

        with open('results.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            headingGrids = ['Time (s)' , 'Left Boundary', ' ']

            gridCounter = 1

            while gridCounter < numHorizontalGridPoints - 1:
                headingGrids.append('')
                headingGrids.append('Gridpoint ' + str(gridCounter))
                headingGrids.append('')
                gridCounter += 1

            headingGrids.append('')
            headingGrids.append('Right Boundary')
            headingGrids.append('')

            writer.writerow(headingGrids)

            #write the next heading row
            gridCounter = 0

            headingTypes = ['']
            while gridCounter < numHorizontalGridPoints:
                headingTypes.append('Depth (m)')
                headingTypes.append('Celerity(m/s)')
                headingTypes.append('Velocity (m/s)')

                gridCounter += 1
            writer.writerow(headingTypes)

            timestepCounter = 0

            for timestep in hartTreeGrid:
                timeValue = timestepCounter * verticalGridSize
                if compareAnalytical:
                    if timeValue % 30 != 0:
                        #only print out 30s timesteps for ease of comparing with analytical solution
                        timestepCounter += 1
                        continue

                gridCounter = 0
                currentTimestepValues = [timeValue]
                for gridpoint in hartTreeGrid[timestepCounter]:
                    currentTimestepValues.append(gridpoint.depth)
                    currentTimestepValues.append(gridpoint.celerity)
                    currentTimestepValues.append(gridpoint.velocity)

                    if diagnosticPrint:
                        print('Grid number ' + str(gridCounter) + ' Celerity is ' + str(gridpoint.celerity)
                              + ' Depth is ' + str(gridpoint.depth) + ' Velocity is ' + str(gridpoint.velocity))

                    gridCounter += 1
                writer.writerow(currentTimestepValues)
                timestepCounter += 1
    else:
        populateDepthBoundary()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
