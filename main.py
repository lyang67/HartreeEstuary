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

def calculateRightBoundaryConditions(time):
    #TODO: implement for real
    rightGridpoint = GridPoint(0,0,0)
    return rightGridpoint

def calculateLeftBoundaryConditions(time):
    #TODO: implement for real
    leftGridpoint = GridPoint(0,0,0)
    return leftGridpoint

def calculateConditionsAtCurrentGridpoint(timestepNumber, gridpointNumber, currentGrid):
    #TODO: implenet this
    print('lol haha')

def populateGrid(mainGrid, timestepSize, xDistanceSize):
    # calculate the number of timesteps we have in a 24 hour period
    numTimesteps = (24 *3600)/timestepSize

    # we have 5 known stations, and 2 boundary conditions, so 5 + 2 grids points in the horizontal direction
    numHorizontalGridPoints = 7

    currentTimestep = 1
    currentGridPoint = 2

    while (currentTimestep <= numTimesteps):
        #do stuff
        leftBoundaryConditions = calculateLeftBoundaryConditions(timestepSize * currentTimestep)
        rightBoundaryConditions = calculateRightBoundaryConditions(timestepSize * currentTimestep)

        #calculate conditions at each grid point between boundary
        #gridpoints 1 and 7 are the boundary conditions points
        while (currentGridPoint < numHorizontalGridPoints):
            #calculate flow conditions at this grid point


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #TODO: the sites don't look 250m apart, but figure this out later
    gridSizeX = 250
    #every 1/2 hour for now, so we don't have to wait forever for the computation to finish while testing
    gridSizeY = 1800
    hartTreeGrid = set_init_conditions()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
