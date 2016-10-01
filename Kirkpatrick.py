import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import collections  as mc

from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from scipy.spatial import delaunay_plot_2d
from scipy.spatial import convex_hull_plot_2d
from random import randint
import numpy as np
import time
import itertools
import Tools as tools

from MyGraph import MyGraph, Piece

# class that performs Kirkpatrick's location method
class Kirkpatrick:
	# max degree allowed to create independent set
	MAX_DEGREE = 8
	# first 3 points in list are the bounding triangle
	POINT_START = 3

	def __init__(self, points):
		# all initial instance variables
		self.points = [tools.roundPoint(point) for point in points]
		# an instance of <MyGraph> used to build DAG location structure
		self.g = None
		# root node of DAG
		self.root = None
		# untouched version of fine graph because <self.g> gets destroyed during building process
		self.untouchedG = None
		# total number of points
		self.N = 0

		# get convex hull of graph before adding bounding triangle
		ch = ConvexHull(self.points)
		interior = [i + Kirkpatrick.POINT_START for i in ch.vertices]

		# add bounding triangle to graph
		boundingTriangle = self.getBoundingTriangle(self.points)
		for point in boundingTriangle:
			self.points.insert(0, tools.roundPoint(point))
	
		# pieces
		pieces = []

		# get exterior triangle
		exterior = [0, 1, 2]

		# get non-interior triangulation and add to pieces collection
		exteriorTri = tools.triangulateRing(self.points, exterior, interior)
		for triangle in exteriorTri:
			pieces.append(Piece(triangle, isLeaf=True, isInside=False))

		# now get interior triangulation and add to pieces collection
		interiorTri = Delaunay(points)
		interiorTri = [[i + Kirkpatrick.POINT_START for i in tri] for tri in interiorTri.simplices]
		for triangle in interiorTri:
			pieces.append(Piece(triangle, isLeaf=True, isInside=True))
		
		# number of points in graph
		self.N = len(self.points)

		# have all the points, now create graph
		self.g = MyGraph(self.points, pieces)
		self.untouchedG = MyGraph(self.points, pieces)

		#self.g.drawMe()
		## keep building DAG until only bounding triangle points are left
		while (self.g.currentN > 3):
			self.root = self.getNextLayer()
			#self.g.drawMe()
		assert (len(self.root) == 1)

		self.root = self.root[0]

	# get bounding triangle to all points on graph and return
	# assumes points do not all lie on a line...
	def getBoundingTriangle(self, pointSet):
		xCoordinates = [row[0] for row in pointSet]
		yCoordinates = [row[1] for row in pointSet]
		
		# get a bounding box
		xMin = float(min(xCoordinates) - 1)
		xMax = float(max(xCoordinates) + 1)
		yMin = float(min(yCoordinates) - 1)
		yMax = float(max(yCoordinates) + 1)

		height = yMax - yMin
		width = xMax - xMin

		# get triangle tip point
		tip = [(xMin + xMax)/2, height + yMax]

		yMin = yMin - height

		# left corner
		boxCorner = [xMin, yMax]
		m = (boxCorner[1] - tip[1])/(boxCorner[0] - tip[0])
		b = boxCorner[1] - m * boxCorner[0]
		leftCorner = [((yMin - b)/m)*1.5, yMin]

		# right corner
		boxCorner = [xMax, yMax]
		m = (boxCorner[1] - tip[1])/(boxCorner[0] - tip[0])
		b = boxCorner[1] - m * boxCorner[0]
		rightCorner = [((yMin - b)/m)*1.5, yMin]

		return [leftCorner, tip, rightCorner]

	# find an independent set using <self.g> and return
	# used by constructor for building kirkpatrick's DAG datastructure
	def findIndependentSet(self):
		# all nodes initially unmarked
		# except for points on bounding triangle
		marked = np.zeros(self.N, dtype = bool)
		for i in xrange(0, Kirkpatrick.POINT_START):
			marked[i] = True

		# mark all nodes with degree greater than <MAX_DEGREE>
		for i in xrange(Kirkpatrick.POINT_START, self.N):
			if not self.g.isActive(i):
				marked[i] = True
			elif self.g.degree(i) > Kirkpatrick.MAX_DEGREE:
				marked[i] = True

		independentSet = []

		# add nodes to independent set and mark it + neighbors
		# keep adding nodes to set until all nodes are marked
		for i in xrange(Kirkpatrick.POINT_START, self.N):
			# if not marked yet...
			if not marked[i]:
				# add to set
				independentSet.append(i)
				# mark it and neighbors
				marked[i] = True
				for neighbor in self.g.getNeighbors(i):
					marked[neighbor] = True

		return independentSet

	# generate the next layer of coarser triangles
	# should not be used if only bounding triangle points are left
	def getNextLayer(self):
		# get the independent set
		indepSet = self.findIndependentSet()

		# for each vertex in indep set..
		for vertex in indepSet:
			# get surrounding polygon
			polygonHole = self.g.getSurroundingPolygon(vertex)

			# triangulate hole and get triangles
			triangulation = tools.triangulatePolygon(self.points, polygonHole)
			
			# create pieces out of triangles
			pieces = []
			for triangle in triangulation:
				pieces.append(Piece(triangle, self.g.getIntersectingPiecesAtP(vertex, triangle)))

			# remove vertex from graph
			self.g.removeVertex(vertex)

			# add new faces to graph
			for piece in pieces:
				self.g.addPiece(piece)

		return pieces

	def drawGraph(self, edges):
		lines = []
		for edge in edges:
			lines.append([self.points[edge[0]], self.points[edge[1]]])
		lc = mc.LineCollection(lines, linewidths=1)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		ax.autoscale()
		ax.margins(0.1)
		plt.show()

	# find triangle point <q> is in
	# if not in bounding triangle, returns a dict with key 'inside' set to False
	# if inside bounding triangle, returns coordinates of triangle in a dict
	# dict also has key 'inside' set to true or false depending on whether
	# located triangle is inside the convex hull of the supplied points
	def locate(self, q):
		traveler = self.root

		# not inside triangle
		if (not tools.insideTriangle(self.points[traveler.p1], self.points[traveler.p2], self.points[traveler.p3], q)):
			return {'inside': False}

		# inside bounding triangle
		# find triangle
		while (not traveler.leaf()):
			found = False
			for child in traveler.children:
				if (tools.insideTriangle(self.points[child.p1], self.points[child.p2], self.points[child.p3], q)):
					traveler = child;
					found = True
					break
			assert (found)

		location = {}
		location['inside'] = traveler.inside()
		location['p1'] = self.points[traveler.p1]
		location['p2'] = self.points[traveler.p2]
		location['p3'] = self.points[traveler.p3]

		return location

	# same functionality as locate() function
	# however animates the query as triangles transition from coarse to fine
	def animatedLocation(self, q):
		traveler = self.root
		# not inside triangle
		if (not tools.insideTriangle(self.points[traveler.p1], self.points[traveler.p2], self.points[traveler.p3], q)):
			return {'inside': False}
			
		# inside bounding triangle
		# find triangle
		while (not traveler.leaf()):
			# create new graph
			fig = plt.figure()
			ax = plt.subplot(111)
			# always draw bounding triangle first
			vertices = [self.points[0], self.points[1], self.points[2], self.points[0]]
			bol = patches.Polygon(vertices, True, fill=False)
			ax.add_patch(bol)
			found = False

			# find child triangle <q> is in
			for child in traveler.children:
				# draw triangle containing <q> as filled
				if (tools.insideTriangle(self.points[child.p1], self.points[child.p2], self.points[child.p3], q)):
					traveler = child;
					found = True
					self._drawPiece(ax, child, True)
				# draw other child triangles as unfilled
				else:
					i = 0
					self._drawPiece(ax, child, False)
			assert (found)

			# finally, draw point
			ax.plot([q[0]], [q[1]], marker='o', color='k', markersize = 3)
			ax.autoscale()
			plt.show()
			assert (found)
		
		location = {}
		location['inside'] = traveler.inside()
		location['p1'] = self.points[traveler.p1]
		location['p2'] = self.points[traveler.p2]
		location['p3'] = self.points[traveler.p3]
		return location

	# leave up to caller to show
	def _drawPiece(self, ax, piece, filled):
		vertices = [self.points[piece.p1], self.points[piece.p2], self.points[piece.p3], self.points[piece.p1]]
		bol = None
		if filled:
			bol = patches.Polygon(vertices, True, fill=True, fc = 'm', ec = 'k')
		else:
			bol = patches.Polygon(vertices, True, fill=False)
		ax.add_patch(bol)

	# draw point on fine graph to show location
	def showPointOnGraph(self, q):
		self.untouchedG.drawMeWithPoint(q)
