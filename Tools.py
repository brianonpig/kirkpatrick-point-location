import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from scipy.spatial import Delaunay
import math
import sys
sys.path.append('./tri-0.2/src/tri/')
from delaunay import ToPointsAndSegments, triangulate
from delaunay import output_triangles, TriangleIterator, InteriorTriangleIterator

# testing method used to translate string of points on a polygon
# into an actual list of points
def processPolygon(string):
	points = []
	trimmed = string[8:-1]
	coordinates = trimmed.split(',')
	for point in coordinates:
		point = point.replace('px', '')
		point = point.split(' ')
		points.append([float(point[0]), float(point[1])])
	return points

# convert a list of items (a, b, c) to a list of pairs ((a,b), (b, c), (c, a))
def listToPairs(items):
	pairs = []
	n = len(items)
	for i in xrange(0, n):
		if i < n - 1:
			pairs.append([items[i], items[i + 1]])
		else:
			pairs.append([items[i], items[0]])
	return pairs

# print a list of points
def printPoints(points):
	i = 0
	for point in points:
		print str(i) + ': ',
		print point
		i+=1

# mod function
def mod(a, b):
	rem = a % b;
	if (rem < 0):
		rem = rem + b;
	return rem;

# round the x and y coordinates of a point
def roundPoint(p):
	return [thouTrunc(p[0]), thouTrunc(p[1])]

# truncate to the thousandths
def thouTrunc(num):
	return math.trunc(num * 1000)/float(1000)

# is the line formed by points <a>, <b>, <c> ccw?
# almost same function as <area>
def ccw(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) > (b[1] - a[1]) * (c[0] - a[0])

# check if two doubles are equal
def doublesEqual(a, b):
	if math.abs(a - b) <= 0.0001:
		return True

# calculates the triangle's size (formed by the "anchor" segment and additional point)
def area2(a, b, c):
     return (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])

# points <a> and <b> forms the anchored segment.
# point <c> is the evaluated point
def isOnLeft(a, b, c):
     return area2(a, b, c) > 0

# points <a> and <b> forms the anchored segment.
# point <c> is the evaluated point
def isOnRight(a, b, c):
     return area2(a, b, c) < 0

# are points <a>, <b>, <c> collinear?
def isCollinear(a, b, c):
	EPSILON = 0.0001
	area = area2(a, b, c)
	if (area >= 0):
		return area <= EPSILON
	return area >= EPSILON

# calculate angle between points <a>, <b>, <c>
def angle(a, b, c):
	val1 = (b[0] - a[0])**2 + (b[1] - a[1])**2;
	val2 = (b[0] - c[0])**2 + (b[1] - c[1])**2;
	val3 = (c[0] - a[0])**2 + (c[1] - a[1])**2;
	return np.arccos( (val1 + val2 - val3) / math.sqrt(4 * val1 * val2) ) * 180/np.pi;

# draw a graph based on points and edges
def drawGraph(points, edges):
	lines = []
	for edge in edges:
		lines.append([points[edge[0]], points[edge[1]]])
	lc = mc.LineCollection(lines, linewidths=1)
	fig, ax = plt.subplots()
	ax.add_collection(lc)
	ax.autoscale()
	ax.margins(0.1)
	plt.show()

# check if <q> lies inside triangle formed by the three points <a>, <b>, <c>
def insideTriangle(a, b, c, q):
	# check if <q> lies on the same side of all 3 lines
	side1 = ccw(a, b, q);
	side2 = ccw(b, c, q);
	side3 = ccw(c, a, q);

	if ((side1 == side2) and (side2 == side3)):
		return True;
	return False;

# do the two line segments <a> and <b> intersect?
def segmentIntersect(a, b):
	assert (len(a) == 2)
	assert (len(b) == 2)

	# check if collinear
	colinear1 = isCollinear(a[0], a[1], b[0])
	colinear2 = isCollinear(a[0], a[1], b[1])
	if (colinear1 and colinear2):
		# get max and min x value
		maxX = max([a[0][0], a[1][0]])
		minX = min([a[0][0], a[1][0]])

		if (b[0][0] >= minX and b[0][0] <= maxX):	
			return True
		elif (b[1][0] >= minX and b[1][0] <= maxX):	
			return True
		return False

	# use line <a> as anchor
	onLeftA0 = isOnLeft(a[0], a[1], b[0])
	onLeftA1 = isOnLeft(a[0], a[1], b[1])
	onRightA0 = isOnRight(a[0], a[1], b[0])
	onRightA1 = isOnRight(a[0], a[1], b[1])
	onePointOnLeftA = onLeftA0 or onLeftA1
	onePointOnRightA = onRightA0 or onRightA1

	if (not (onePointOnLeftA and onePointOnRightA)):
		return False

	# use other line <b> as anchor
	onLeftB0 = isOnLeft(b[0], b[1], a[0])
	onLeftB1 = isOnLeft(b[0], b[1], a[1])
	onRightB0 = isOnRight(b[0], b[1], a[0])
	onRightB1 = isOnRight(b[0], b[1], a[1])
	onePointOnLeftB = onLeftB0 or onLeftB1
	onePointOnRightB = onRightB0 or onRightB1

	if (not (onePointOnLeftB and onePointOnRightB)):
		return False

	return True

# triangulate a polygon with no interior points
# returns triangles as indices into <points>
def triangulatePolygon(points, polygon):
	# each point becomes a tuple
	pointTuples = [(point[0], point[1]) for point in points]
	pointDict = {}
	for i, point in enumerate(pointTuples):
		pointDict[point] = i

	polygonVals = [pointTuples[i] for i in polygon]
	polygonVals.append(polygonVals[0])

	pts_segs = ToPointsAndSegments()
	pts_segs.add_polygon([polygonVals])
	dt = triangulate(pts_segs.points, pts_segs.infos, pts_segs.segments)

	triangles = []

	# retrieve triangles
	for t in InteriorTriangleIterator(dt):
		triangle = []
		for vertex in t.vertices:
			triangle.append(pointDict[(vertex.x, vertex.y)])
		triangles.append(triangle)

	return triangles

# triangulate a space between an exterior polygon <exterior>
# and interior polygon <interior>
# returns triangles as indices into <points>
def triangulateRing(points, exterior, interior):
	# each point becomes a tuple
	pointTuples = [(point[0], point[1]) for point in points]
	pointDict = {}
	for i, point in enumerate(pointTuples):
		pointDict[point] = i

	exteriorVals = [pointTuples[i] for i in exterior]
	exteriorVals.append(pointTuples[exterior[0]])
	interiorVals = [pointTuples[i] for i in interior]
	interiorVals.append(pointTuples[interior[0]])

	# create points and segments for triangulation
	pts_segs = ToPointsAndSegments()
	pts_segs.add_polygon([exteriorVals, interiorVals])

	dt = triangulate(pts_segs.points, pts_segs.infos, pts_segs.segments)

	triangles = []

	# retrieve triangles
	for t in InteriorTriangleIterator(dt):
		triangle = []
		for vertex in t.vertices:
			triangle.append(pointDict[(vertex.x, vertex.y)])
		triangles.append(triangle)

	# graph triangles
	edges = []
	for triangle in triangles:
		for edge in listToPairs(triangle):
			edges.append(edge)

	return triangles

def main():

	# test triangulation function
	points = [(0, 0), (0,5), (5,5), (5, 0), (1, 1), (1, 2), (2, 2), (2, 1)]
	exterior = [0, 1, 2, 3]
	interior = [4, 5, 6, 7]

	print triangulateRing(points, exterior, interior)

if __name__ == "__main__":
    main()