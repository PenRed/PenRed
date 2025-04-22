import bpy
import gpu
from gpu_extras.batch import batch_for_shader
from mathutils import Vector
from mathutils import Color
from math import log

def clamp(v, minV, maxV):
    return max(minV, min(v, maxV))

defaultMaxBatches = 50

class TrackManager:
    def __init__(self, maxBatches=defaultMaxBatches):
        self.batches = []
        self.shader = gpu.shader.from_builtin('SMOOTH_COLOR')
        self.maxBatches = maxBatches  # Adjust based on your needs
        self.energyColors = (Color((1.0,1.0,0.0)), Color((1.0,0.0,0.0)))
        self.timeColors = (Color((0.8,0.8,1.0)), Color((0.0,0.0,8.0)))

        self.energyRange = (50.0,1.0e6)
        self.timeRange = (0.0,30.0)

        self.eLogScale = False
        self.outOfRange = True
        
        self.colorType = "ENERGY"

    def setOutOfRange(self, enable):
        if self.outOfRange == enable:
            return False
        self.outOfRange = enable
        return True

    def setColorType(self, newColorType):
        if self.colorType == newColorType:
            return False
        self.colorType = newColorType
        return True
        
    def setELogScale(self, enable):
        if self.eLogScale == enable:
            return False
        self.eLogScale = enable
        return self.colorType == "ENERGY"

    def setERange(self, newRange):
        if newRange[0] == self.energyRange[0] and newRange[1] == self.energyRange[1]:
            return False
        self.energyRange =  list(newRange)
        self.energyRange[0] = max(self.energyRange[0], 1.0)
        self.energyRange[1] = max(self.energyRange[0], self.energyRange[1])
        return self.colorType == "ENERGY"

    def setTRange(self, newRange):
        if newRange[0] == self.timeRange[0] and newRange[1] == self.timeRange[1]:
            return False
        self.timeRange =  list(newRange)
        self.timeRange[0] = max(self.timeRange[0], 0.0)
        self.timeRange[1] = max(self.timeRange[0], self.timeRange[1])
        return self.colorType == "TIME"
    
    def fitERange(self, newRange):
        updated = False
        if self.energyRange[0] > newRange[0]:
            self.energyRange[0] = newRange[0]
            updated = True
        if self.energyRange[1] < newRange[1]:
            self.energyRange[1] =  newRange[1]
            updated = True

        if updated:
            return self.colorType == "ENERGY"
        return False

    def fitTRange(self, newRange):
        updated = False
        if self.timeRange[0] > newRange[0]:
            self.timeRange[0] = newRange[0]
            updated = True
        if self.timeRange[1] < newRange[1]:
            self.timeRange[1] =  newRange[1]
            updated = True

        if updated:
            return self.colorType == "TIME"
        return False
    
    def draw(self, enabledTracks):
        if len(self.batches) > 0:
            self.shader.bind()
            n = min(len(enabledTracks), len(self.batches))
            i = 0
            while i < n:
                if enabledTracks[i]:
                    self.batches[i].draw(self.shader)
                i = i + 1

    def removeBatch(self, i):
        del self.batches[i]
                
    def addBatch(self, tracks):

        if len(tracks) == 0:
            return
        if len(self.batches) >= self.maxBatches:
            return

        vertices = []
        colors = []
        indices = []
        for track in tracks:
            
            if self.colorType == "ENERGY":
                
                minColor = self.energyColors[0]
                maxColor = self.energyColors[1]
                
                if self.eLogScale:
                    minVal = log(self.energyRange[0])
                    maxVal = log(self.energyRange[1])
                    dVal = maxVal - minVal
                    if dVal == 0:
                        dVal = 1.0
                    colFactors = [(log(point.energy) - minVal)/dVal for point in track.points]
                else:
                    minVal = self.energyRange[0]
                    maxVal = self.energyRange[1]
                    dVal = maxVal - minVal
                    if dVal == 0:
                        dVal = 1.0                    
                    colFactors = [(point.energy - minVal)/dVal for point in track.points]
                    
            elif self.colorType == "TIME":
                minVal = self.timeRange[0]
                maxVal = self.timeRange[1]
                dVal = maxVal - minVal
                if dVal == 0:
                    dVal = 1.0
                
                colFactors = [(point.time - minVal)/dVal for point in track.points]
            else:
                return

            # Check if points out of range have to be drawn
            if self.outOfRange:
                #Clamp color factors
                colFactors = [clamp(f, 0.0, 1.0) for f in colFactors]
                trackVert = [point.position for point in track.points]
            else:
                #Select only points in range
                inRangeIndex = [i for i, x in enumerate(colFactors) if 0.0 <= x <= 1.0]

                colFactors = [colFactors[i] for i in inRangeIndex]
                trackVert = [track.points[i].position for i in inRangeIndex]

            initIndex = len(vertices)
            if len(trackVert) < 2:
                continue
                
            # Extend vertice vector
            vertices.extend(trackVert)

            trackCol = [minColor*(1.0-f) + maxColor*f for f in colFactors]
            colors.extend([(*c, 1.0) for c in trackCol])

            trackIndex = [(initIndex + i, initIndex + i + 1) for i in range(len(trackVert)-1)]
            indices.extend(trackIndex)

        # Create batch
        batch = batch_for_shader(
            self.shader, 'LINES',
            {"pos": vertices, "color": colors},
            indices=indices
        )

        self.batches.append(batch)
        
    def updateTracks(self, trackFiles):

        # Clear old batch tracks
        self.clear()

        # Create batches
        ifile = 0
        limit = min(len(trackFiles), self.maxBatches)
        while ifile < limit:
            self.addBatch(trackFiles[ifile].tracks)
            ifile = ifile + 1
            
    def clear(self):
        self.batches.clear()

# Singleton instance
_trackManagerInstance = None

def getTrackManager(maxBatches=defaultMaxBatches):
    global _trackManagerInstance
    if _trackManagerInstance is None:
        _trackManagerInstance =  TrackManager(maxBatches)
    return _trackManagerInstance

def clearTrackManager():
    global _trackManagerInstance
    if _trackManagerInstance is not None:
        _trackManagerInstance.clear()
        del _trackManagerInstance
        _trackManagerInstance = None
