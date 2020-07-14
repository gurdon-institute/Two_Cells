# Selects best focussed slice and segments 2 cells from C2 (DIC) in a 5D stack. Segments a nucleus in each cell if one is labelled in C1.
# Measures C1 mean intensity over time for each cell and nucleus.
#
# 	- by Richard Butler, Gurdon Institute Imaging Facility


import itertools
import math as maths

from java.awt import Color, Rectangle, BasicStroke, Dimension
from java.awt.geom import Ellipse2D, Path2D

from ij import IJ, ImagePlus, ImageStack
from ij.plugin.filter import MaximumFinder, RankFilters, ThresholdToSelection, EDM, BackgroundSubtracter
from ij.process import ImageStatistics, StackStatistics, Blitter, ImageProcessor, ByteProcessor, FloatProcessor, AutoThresholder, FloodFiller, EllipseFitter
from ij.measure import Measurements, ResultsTable
from ij.gui import Roi, ShapeRoi, PolygonRoi, Line, Overlay

from org.jfree.chart import JFreeChart, ChartFactory, ChartFrame, LegendItemCollection, LegendItem
from org.jfree.chart.plot import PlotOrientation
from org.jfree.data.xy import DefaultXYDataset


def fillHoles(ip):
    width = ip.getWidth()
    height = ip.getHeight()
    ff = FloodFiller(ip)
    ip.setColor(127)
    foreground = 127
    background = 0
    for y in range(height):
        if ip.getPixel(0,y)==background:
            ff.fill(0, y)
        if ip.getPixel(width-1,y)==background:
            ff.fill(width-1, y)
    for x in range(width):
        if ip.getPixel(x,0)==background:
            ff.fill(x, 0)
        if ip.getPixel(x,height-1)==background:
            ff.fill(x, height-1)
    n = width*height
    for i in range(n):
        if ip.get(i)==127:
            ip.set(i, 0)
        else:
            ip.set(i, 255)

def watershed(ip):
    TOL = 0.5
    floatEdm = EDM().makeFloatEDM(ip, 0, False)
    maxIp = MaximumFinder().findMaxima(floatEdm, TOL, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
    if (maxIp != None):
        ip.copyBits(maxIp, 0, 0, Blitter.AND)

def onEdge(roi):
    for xy in range(max(W,H)):
        if roi.contains(0,xy) or roi.contains(xy,0) or roi.contains(W-1,xy) or roi.contains(xy,H-1):
            return True
    return False

def getDICfocus(imp):
    stack = imp.getStack()
    focusStack = ImageStack(W,H)
    for t in range(1,T+1):
        sd = [0 for z in range(Z+1)]
        for z in range(1,Z+1):	#get best focussed C2 (DIC) slices
            ip = stack.getProcessor(imp.getStackIndex(2,z,t)).convertToFloat()
            ip.findEdges()
            sd[z] = ImageStatistics.getStatistics(ip, Measurements.STD_DEV, cal).stdDev
        focusZ = sd.index(max(sd))
        focusSlice = stack.getProcessor(imp.getStackIndex(2,focusZ,t)).convertToFloatProcessor()
        focusStack.addSlice("Z"+str(focusZ), focusSlice)

    #ImagePlus("DIC focus", focusStack).show()
    #exit(0)

    return focusStack

def getCells(dicStack):
    outStack = ImageStack(W,H)

    cells = [None for t in range(T+1)]

    for t in range(1,T+1):
        mapp = dicStack.getProcessor(t).convertToFloatProcessor()

        mapp.subtract( mapp.getStatistics().mean )
        mapp.abs()

        RankFilters().rank(mapp, 1.0, RankFilters.VARIANCE)
        mapp.sqrt()

        mapp.blurGaussian(5)

        hist = mapp.getHistogram(256)
        stats = mapp.getStatistics()

        thresh = AutoThresholder().getThreshold( AutoThresholder.Method.Otsu, hist )
        thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min

        mask = ByteProcessor(W,H)
        for i in range(W*H):
            value = mapp.getf(i)
            bite = 255 if value>=thresh else 0
            mask.set(i, bite)

        fillHoles(mask)
        ed = 3
        for e in range(ed): mask.erode(1, 0)
        for d in range(ed): mask.dilate(1, 0)

        watershed(mask)

        minA = 5000 #px²

        mask.setThreshold(255,255, ImageProcessor.NO_LUT_UPDATE)
        composite = ThresholdToSelection().convert(mask)

        rois = ShapeRoi(composite).getRois()
        keep = []
        for roi in rois:
            if roi.getStatistics().area >= minA:
                if not onEdge(roi):
                    keep.append(roi)
                else:
                    edgeRoi = ShapeRoi(roi)
                    edgeRoi.setPosition(0,0,t)
                    edgeRoi.setStrokeColor(Color.YELLOW)
                    ol.add(edgeRoi)
        print("T"+str(t)+" using "+str(len(keep))+"/"+str(len(rois))+" ROIs")
        rois = keep
        #rois = [ roi for roi in rois if roi.getStatistics().area >= minA and not onEdge(roi) ]	#keep big enough and not on edges

        # if there is only one Roi, cut it along the fitted ellipse minor axis
        if len(rois)==1:
            el = EllipseFitter()
            mask.setRoi(rois[0])
            el.fit(mask, None)
            el.makeRoi(mask)
            theta = el.angle * (maths.pi/180.0)

            length = el.major/2.0
            dy = maths.sin(theta)* length
            dx = maths.cos(theta)* length

            #major axis
            lineX0 = el.xCenter - dx
            lineY0 = el.yCenter + dy
            lineX1 = el.xCenter + dx
            lineY1 = el.yCenter - dy
            line = Line(lineX0, lineY0, lineX1, lineY1)
            line.setStrokeColor(Color.BLUE)
            line.setStrokeWidth(1)
            line.setPosition(0,0,t)
            ol.add(line)

            #minor axis scaled length to make sure cut ends are outside Roi
            cutX0 = el.xCenter + dy*100
            cutY0 = el.xCenter + dx*100
            cutX1 = el.yCenter - dy*100
            cutY1 = el.yCenter - dx*100

            cut = Line(cutX0,cutY0, cutX1, cutY1)
            cut.setStrokeWidth(2)
            cut = PolygonRoi( cut.getFloatPolygon(), PolygonRoi.POLYGON )

            mask.setColor(0)
            mask.fill(cut)
            composite = ThresholdToSelection().convert(mask)

            rois = ShapeRoi(composite).getRois()
            rois = [ roi for roi in rois if roi.getStatistics().area >= minA ]
        print(str(t) + ":" + str(len(rois)))

        rois = [ PolygonRoi(roi.getInterpolatedPolygon(20, True), PolygonRoi.POLYGON) for roi in rois ]
        rois = [ PolygonRoi(roi.getConvexHull(), PolygonRoi.POLYGON) for roi in rois ]

        rois = sorted(list(rois), key=lambda roi:roi.getLength() )	#size order
        rois = rois[-2:]											#keep 2 biggest
        rois = sorted(list(rois), key=lambda roi:roi.getStatistics().xCentroid+roi.getStatistics().yCentroid )	#top left to bottom right order

        if len(rois)>0:
            rois[0].setStrokeColor(Color.RED)
            rois[0].setPosition(0, 0, t)
            ol.add(rois[0])
        if len(rois)>1:
            rois[1].setStrokeColor(Color.GREEN)
            rois[1].setPosition(0, 0, t)
            ol.add(rois[1])
            cells[t] = (rois[0], rois[1])


    return cells

def getC1Projection(imp, norm):
    stack = imp.getStack()
    projStack = ImageStack(W,H)
    for t in range(1,T+1):
        proj = FloatProcessor(W,H)
        for z in range(1,Z+1):	#max intensity projection
            ip = stack.getProcessor(imp.getStackIndex(1,z,t)).convertToFloatProcessor()
            proj.copyBits(ip, 0,0, Blitter.MAX)
        projStack.addSlice(proj)

    if norm:	#Z normalise
        stats = StackStatistics(ImagePlus("wrapper",projStack))
        for t in range(1,T+1):
            ip = projStack.getProcessor(t)
            ip.subtract(stats.mean)
            ip.multiply( 1.0/(stats.stdDev) )

    return projStack

def getNuclei(stack):

    nuclei = [ [] for t in range(T+1) ]

    minNucleusA = 50.0			#µm²
    maxNucleusA = 500.0
    sigma = 0.5 * maths.sqrt(minNucleusA/maths.pi) / cal.pixelWidth	#px
    k = 5

    for t in range(1,T+1):
        proc = stack.getProcessor(t).duplicate()
        sub = proc.duplicate()

        proc.blurGaussian(sigma)
        sub.blurGaussian(sigma*k)
        proc.copyBits(sub, 0,0, Blitter.SUBTRACT)
        hist = proc.getHistogram(256)
        stats = proc.getStatistics()
        thresh = AutoThresholder().getThreshold( AutoThresholder.Method.MaxEntropy, hist )
        thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
        proc.setThreshold(thresh, 99999999, ImageProcessor.NO_LUT_UPDATE)
        composite = ThresholdToSelection().convert(proc)
        rois = ShapeRoi(composite).getRois()
        for nuc in rois:
            proc.setRoi(nuc)
            if proc.getStatistics().mean < thresh: continue	#exclude composite holes
            area = nuc.getStatistics().area * cal.pixelWidth * cal.pixelHeight
            if area >= minNucleusA and area <= maxNucleusA:
                circ = 4*maths.pi*(area/pow(nuc.getLength()*cal.pixelWidth, 2))
                if circ >= 0.65:
                    nuclei[t].append(nuc)

        nuclei[t] = sorted( list(nuclei[t]), key=lambda nuc:nuc.getLength(), reverse=True ) #largest to smallest

    return nuclei

def measure(stack, cells, nuclei):
    time = [ (t-1)*cal.frameInterval for t in range(T+1) ]
    cellValues0 = [ 0.0 for t in range(T+1) ]
    cellValues1 = [ 0.0 for t in range(T+1) ]
    cellAreas0 = [ 0.0 for t in range(T+1) ]
    cellAreas1 = [ 0.0 for t in range(T+1) ]
    nucleusValues0 = [ 0.0 for t in range(T+1) ]
    nucleusValues1 = [ 0.0 for t in range(T+1) ]
    nucleusAreas0 = [ 0.0 for t in range(T+1) ]
    nucleusAreas1 = [ 0.0 for t in range(T+1) ]
    nonNucleusValues0 = [ 0.0 for t in range(T+1) ]
    nonNucleusValues1 = [ 0.0 for t in range(T+1) ]

    for t in range(1,T+1):
        ip = stack.getProcessor(t)

        if cells[t] is None:
            continue


        #subtract background Z from all intensity Z measurements
        if cells [t] is None:
            print("Nocellsfound" + str(t))
        bothCells = ShapeRoi(cells[t][0]).or(ShapeRoi(cells[t][1]))
        backRoi = ShapeRoi(Rectangle(0,0,imp.getWidth(),imp.getHeight())).not( bothCells )


        ip.setRoi(backRoi)
        backMean = ip.getStatistics().mean

        ip.setRoi( cells[t][0] )
        stats0 = ip.getStatistics()
        cellValues0[t] = stats0.mean - backMean
        cellAreas0[t] = stats0.area * cal.pixelWidth * cal.pixelHeight
        nuc0 = None
        for nuc in nuclei[t]:
            rect = nuc.getBounds()
            nx = int(rect.x+(rect.width/2.0))
            ny = int(rect.y+(rect.height/2.0))
            if cells[t][0].contains(nx,ny):
                nuc0 = nuc
                break
        if nuc0 is not None:
            ip.setRoi( nuc0 )
            nucStats0 = ip.getStatistics()
            nucleusValues0[t] = nucStats0.mean - backMean
            nucleusAreas0[t] = nucStats0.area * cal.pixelWidth * cal.pixelHeight
            nuc0.setPosition(0,0,t)
            nuc0.setStrokeColor(Color.CYAN)
            ol.add(nuc0)
            nonnucRoi0 = ShapeRoi(cells[t][0]).not( ShapeRoi(nuc0) )
            ip.setRoi( nonnucRoi0 )
            nonNucleusValues0[t] = ip.getStatistics().mean - backMean

        ip.setRoi( cells[t][1] )
        stats1 = ip.getStatistics()
        cellValues1[t] = stats1.mean - backMean
        cellAreas1[t] = stats1.area * cal.pixelWidth * cal.pixelHeight
        nuc1 = None
        for nuc in nuclei[t]:
            rect = nuc.getBounds()
            nx = int(rect.x+(rect.width/2.0))
            ny = int(rect.y+(rect.height/2.0))
            if cells[t][1].contains(nx,ny):
                nuc1 = nuc
                break
        if nuc1 is not None:
            ip.setRoi( nuc1 )
            nucStats1 = ip.getStatistics()
            nucleusValues1[t] = nucStats1.mean - backMean
            nucleusAreas1[t] = nucStats1.area * cal.pixelWidth * cal.pixelHeight
            nuc1.setPosition(0,0,t)
            nuc1.setStrokeColor(Color.CYAN)
            ol.add(nuc1)
            nonnucRoi1 = ShapeRoi(cells[t][1]).not( ShapeRoi(nuc1) )
            ip.setRoi( nonnucRoi1 )
            nonNucleusValues1[t] = ip.getStatistics().mean - backMean

    rt = ResultsTable()
    rt.showRowNumbers(False)
    for t in range(1,T+1):
        rt.setValue("Time ("+cal.getTimeUnit()+")", t-1, IJ.d2s(time[t],1))
        areaRatio = cellAreas0[t] / cellAreas1[t] if cellAreas0[t]>0 and cellAreas1[t]>0 else 0.0
        rt.setValue("Cell 0:Cell 1 Area Ratio", t-1, areaRatio)

        nucleusRatio = nucleusValues0[t] / nucleusValues1[t] if nucleusValues0[t]>0 and nucleusValues1[t]>0 else 0.0
        rt.setValue("Cell 0:Cell 1 Nucleus Ratio", t-1, nucleusRatio)
        nonNucleusRatio = nonNucleusValues0[t] / nonNucleusValues1[t] if nonNucleusValues0[t]>0 and nonNucleusValues1[t]>0 else 0.0
        rt.setValue("Cell 0:Cell 1 Non-Nucleus Ratio", t-1, nonNucleusRatio)

        nnnRatio0 = nucleusValues0[t] / nonNucleusValues0[t] if nucleusValues0[t]>0 and nonNucleusValues0[t]>0 else 0.0
        rt.setValue("Cell 0 Nucleus:Non-Nucleus Ratio", t-1, nnnRatio0)
        nnnRatio1 = nucleusValues1[t] / nonNucleusValues1[t] if nucleusValues1[t]>0 and nonNucleusValues1[t]>0 else 0.0
        rt.setValue("Cell 1 Nucleus:Non-Nucleus Ratio", t-1, nnnRatio1)

        rt.setValue("Cell 0 (red) Area ("+cal.getUnit()+u"\u00b2"+")", t-1, cellAreas0[t])
        rt.setValue("Cell 0 Nucleus Area ("+cal.getUnit()+u"\u00b2"+")", t-1, nucleusAreas0[t])
        rt.setValue("Cell 0 All", t-1, cellValues0[t])
        rt.setValue("Cell 0 Nucleus", t-1, nucleusValues0[t])
        rt.setValue("Cell 0 Non-Nucleus", t-1, nonNucleusValues0[t])
        rt.setValue("Cell 1 (green) Area ("+cal.getUnit()+u"\u00b2"+")", t-1, cellAreas1[t])
        rt.setValue("Cell 1 Nucleus Area ("+cal.getUnit()+u"\u00b2"+")", t-1, nucleusAreas1[t])
        rt.setValue("Cell 1 All", t-1, cellValues1[t])
        rt.setValue("Cell 1 Nucleus", t-1, nucleusValues1[t])
        rt.setValue("Cell 1 Non-Nucleus", t-1, nonNucleusValues1[t])
    rt.show(imp.getTitle()+"-Results")

    dataset = DefaultXYDataset()
    dataset.addSeries( "Cell 0", [time[1:], cellValues0[1:]] )
    dataset.addSeries( "Cell 1", [time[1:], cellValues1[1:]] )
    dataset.addSeries( "Nucleus 0", [time[1:], nucleusValues0[1:]] )
    dataset.addSeries( "Nucleus 1", [time[1:], nucleusValues1[1:]] )
    dataset.addSeries( "Non-Nucleus 0", [time[1:], nonNucleusValues0[1:]] )
    dataset.addSeries( "Non-Nucleus 1", [time[1:], nonNucleusValues1[1:]] )

    chart = ChartFactory.createScatterPlot( imp.getTitle(), "Time ("+cal.getTimeUnit()+")", "Intensity Z", dataset, PlotOrientation.VERTICAL, True,True,False )
    plot = chart.getPlot()

    plot.setBackgroundPaint(Color(64, 128, 255))
    plot.setDomainGridlinePaint(Color.BLACK)
    plot.setRangeGridlinePaint(Color.BLACK)

    renderer = plot.getRenderer()
    legend = LegendItemCollection()
    shapeR = 2.0
    nucShape = Ellipse2D.Float(-shapeR,-shapeR,shapeR*2,shapeR*2)
    nonNucShape = Path2D.Float()
    nonNucShape.moveTo(-shapeR,-shapeR)
    nonNucShape.lineTo(shapeR,shapeR)
    nonNucShape.moveTo(shapeR,-shapeR)
    nonNucShape.lineTo(-shapeR,shapeR)
    for s in range(dataset.getSeriesCount()):

        if s == 0:
            renderer.setSeriesLinesVisible(s, True)
            renderer.setSeriesShapesVisible(s, False)
            renderer.setSeriesStroke(s, BasicStroke(1))
            renderer.setSeriesPaint(s, Color.RED)
            legend.add( LegendItem("Cell 0", Color.RED) )
        elif s == 1:
            renderer.setSeriesLinesVisible(s, True)
            renderer.setSeriesShapesVisible(s, False)
            renderer.setSeriesStroke(s, BasicStroke(1))
            renderer.setSeriesPaint(s, Color.GREEN)
            legend.add( LegendItem("Cell 1", Color.GREEN) )
        elif s == 2:
            renderer.setSeriesLinesVisible(s, False)
            renderer.setSeriesShapesVisible(s, True)
            renderer.setSeriesShape(s, nucShape)
            renderer.setSeriesPaint(s, Color.RED)
        elif s == 3:
            renderer.setSeriesLinesVisible(s, False)
            renderer.setSeriesShapesVisible(s, True)
            renderer.setSeriesShape(s, nucShape)
            renderer.setSeriesPaint(s, Color.GREEN)
        elif s == 4:
            renderer.setSeriesLinesVisible(s, False)
            renderer.setSeriesShapesVisible(s, True)
            renderer.setSeriesShape(s, nonNucShape)
            renderer.setSeriesPaint(s, Color.RED)
        elif s == 5:
            renderer.setSeriesLinesVisible(s, False)
            renderer.setSeriesShapesVisible(s, True)
            renderer.setSeriesShape(s, nonNucShape)
            renderer.setSeriesPaint(s, Color.GREEN)


    plot.setFixedLegendItems(legend)

    frame = ChartFrame(imp.getTitle()+" Z-Normalised Intensity", chart)
    frame.pack()
    frame.setSize( Dimension(800, 800) )
    frame.setLocationRelativeTo(None)
    frame.setVisible(True)


imp = IJ.getImage()
cal = imp.getCalibration()
if cal.pixelWidth > 0.5:
    IJ.error("Bad calibration: "+str(cal.pixelWidth)+" "+cal.getUnit())
    exit(0)
ol = Overlay()
stack = imp.getStack()
W = imp.getWidth()
H = imp.getHeight()
Z = imp.getNSlices()
T = imp.getNFrames()


focusStack = getDICfocus(imp)
projC1 = getC1Projection(imp, True)

cells = getCells(focusStack)

nuclei = getNuclei(projC1)

measure(projC1, cells, nuclei)
imp.setOverlay(ol)
