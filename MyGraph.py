import numpy as np
import matplotlib.pyplot as plt
import Tools as tools
from matplotlib import collections  as mc

# class used to represent a triangle (with pointers) in the Kirkpatrick point location algorithm
# a Piece must be associated with a fixed list of points to make sense, as its triangle
# is defined by indices into that list
class Piece:
	def __init__(self, triangle, children=None, isLeaf=False, isInside=False):
		assert (len(triangle) == 3)

		triangle.sort()

		# indices representing triangle points
		self.p1 = triangle[0]
		self.p2 = triangle[1]
		self.p3 = triangle[2]
		# list of children (intersecting triangles in finer level)
		self.children = children
		# is this piece a leaf?
		self.isLeaf = isLeaf
		# is this piece inside the original convex hull of the points?
		self.isInside = isInside 

	# set the children of this piece
	def setChildren(children):
		self.children = children

	# does this piece represent the same triangle as <otherPiece>?
	def equals(self, otherPiece):
		if (self.p1 != otherPiece.p1):
			return False
		elif (self.p2 != otherPiece.p2):
			return False
		elif (self.p3 != otherPiece.p3):
			return False
		return True

	# is this piece a leaf?
	def leaf(self):
		return self.isLeaf

	# is this piece inside the original convex hull of the points?
	def inside(self):
		assert(self.leaf())
		return self.isInside

	# string representation of piece (aka the triangle indices)
	def toString(self):
		string = "piece: " + str(self.p1) +", " + str(self.p2) +", "+ str(self.p3)
		return string

# class used to represent the graph used to generate the datastructure for 
# Kirkpatrick location. It essentially uses an adjacency matrix representation
# for a graph. Instead of simply storing a boolean to determine whether an edge exists, however,
# it stores the Piece (the triangle(s)) that that edge is a part of. Consequently, the 
# "adjacency matrix" is an N x N x 2 graph (2 faces possible at an edge), where N = # total points. 
class MyGraph:

	# number points remain constant
	def __init__(self, points, pieces):
		# all possible points
		self.points = points
		# number of max number of points
		self.N = len(points)
		# number of current number of points
		self.currentN = len(points)
		# a collection of triangles of the current layer
		self.pieces = pieces
		# an array of which points are active (not removed yet)
		self.activePoints = np.ones(self.N, dtype = bool)
		# an array of pieces ordered by edge
		self.edgeToFace = np.empty((self.N, self.N, 2), dtype=object)
		self.edgeToFace.fill(None)

		# add pieces
		for piece in self.pieces:
			self.addPiece(piece)

	# add a piece to this graph
	def addPiece(self, piece):

		# three points
		a = piece.p1
		b = piece.p2
		c = piece.p3

		# add face to each edge
		self._addFaceToEdge(a, b, piece)
		self._addFaceToEdge(a, c, piece)
		self._addFaceToEdge(b, c, piece)

	# erase vertex and associated faces
	def removeVertex(self, vertex):
		assert (vertex >= 3) and (vertex < self.N)
		assert (self.activePoints[vertex])

		faces = self.getFacesAtVertex(vertex)
		for face in faces:
			self.removeFace(face)

		self.activePoints[vertex] = False
		self.currentN -=1

	# remove piece (face) <piece> from this graph
	def removeFace(self, piece):
		self._removeFaceAtEdge(piece, [piece.p1, piece.p2])
		self._removeFaceAtEdge(piece, [piece.p2, piece.p3])
		self._removeFaceAtEdge(piece, [piece.p3, piece.p1])

	# remove a piece <piece> at a specific edge <edge> if it exists
	def _removeFaceAtEdge(self, piece, edge):
		# do it one way first
		twinFaces = self.edgeToFace[edge[0]][edge[1]]
		for i, poss in enumerate(twinFaces):
			if not (poss is None):
				if poss.equals(piece):
					self.edgeToFace[edge[0]][edge[1]][i] = None

		# do it other way now
		twinFaces = self.edgeToFace[edge[1]][edge[0]]
		eliminated = False
		for i, poss in enumerate(twinFaces):
			if not (poss is None):
				if poss.equals(piece):
					self.edgeToFace[edge[1]][edge[0]][i] = None

	# what is the degree of <vertex> in graph?
	# cannot ask for degree of deleted vertex
	def degree(self, vertex):
		assert (vertex >= 3) and (vertex < self.N)
		assert (self.activePoints[vertex])

		count = 0
		for i in xrange(0, self.N):
			if self.edgeToFace[vertex][i][0] is not None:
				count +=1
		return count

	# has <vertex> been removed yet?
	def isActive(self, vertex):
		return self.activePoints[vertex]

	# get neighbors of <vertex>
	def getNeighbors(self, vertex):
		assert (vertex >= 3) and (vertex < self.N)
		assert (self.activePoints[vertex])
		neighbors = []
		for i in xrange(0, self.N):
			if self.edgeToFace[vertex][i][0] is not None:
				neighbors.append(i)
			elif self.edgeToFace[vertex][i][1] is not None:
				neighbors.append(i)
		return neighbors

	# get the (at most 2) pieces at edge <edge>
	def getPieces(self, edge):
		assert (len(edge) == 2)
		return self.edgeToFace[edge[0]][edge[1]]

	# given an active point <vertex>, return the polygon "hole" that would exist
	# if <vertex> were to be deleted. The points are returned in counterclockwise order.
	def getSurroundingPolygon(self, vertex):
		assert (vertex >= 3) and (vertex < self.N)
		assert (self.activePoints[vertex])

		polygon = []

		# get neighbors
		neighbors = self.getNeighbors(vertex)

		# use the first neighbor as a starting point
		edge = [vertex, neighbors[0]]

		# get starting and ending faces
		faces = self.getPieces(edge)

		traveler = self._getCounterClockwiseFace(faces, edge)
		endingFace = self._getOtherFace(faces, traveler)

		polygon.append(edge[1])

		# get all other faces in order
		while (traveler != endingFace):
			edge[1] = self._getOtherPoint(traveler, vertex, edge[1])
			traveler = self._getOtherFace(self.getPieces(edge), traveler)
			polygon.append(edge[1])

		return polygon

	# get the pieces at vertex <vertex>
	def getFacesAtVertex(self, vertex):
		allFaces = []
		# get all neighbors of vertex
		neighbors = self.getNeighbors(vertex)

		for neighbor in neighbors:
			allFaces.append(self.edgeToFace[vertex][neighbor][0])
			allFaces.append(self.edgeToFace[vertex][neighbor][1])
		return set(allFaces)

	# find all pieces at <vertex> that intersect <triangle>
	def getIntersectingPiecesAtP(self, vertex, triangle):
		intersections = []

		l1 = [self.points[triangle[0]], self.points[triangle[1]]]
		l2 = [self.points[triangle[1]], self.points[triangle[2]]]
		l3 = [self.points[triangle[2]], self.points[triangle[0]]]
		segments = [l1, l2, l3]

		# get all neighbors of vertex
		neighbors = self.getNeighbors(vertex)

		# if point is in polygon it intersects all old triangles
		if (tools.insideTriangle(self.points[triangle[0]], self.points[triangle[1]], 
		self.points[triangle[2]], self.points[vertex])):
			return self.getFacesAtVertex(vertex)

		# for each edge... check intersection with any legs of triangle
		for neighbor in neighbors:
			edge = [self.points[vertex], self.points[neighbor]]
			for segment in segments:
				# if they intersect
				# add appropriate faces
				if (tools.segmentIntersect(edge, segment)):
					intersections.append(self.edgeToFace[vertex][neighbor][0])
					intersections.append(self.edgeToFace[vertex][neighbor][1])

		return set(intersections)

	# get current number of points
	def getCurrentNumPoints(self):
		return self.currentN

	# add <piece> to edge
	def _addFaceToEdge(self, a, b, piece):
		# one way edge
		if self.edgeToFace[a][b][0] is None:
			self.edgeToFace[a][b][0] = piece
		else:
			self.edgeToFace[a][b][1] = piece

		# add for other way
		if self.edgeToFace[b][a][0] is None:
			self.edgeToFace[b][a][0] = piece
		else:
			self.edgeToFace[b][a][1] = piece

	# return other point on a piece
	def _getOtherPoint(self, piece, a, b):
		if (piece.p1 != a) and (piece.p1 != b):
			return piece.p1
		elif (piece.p2 != a) and (piece.p2 != b):
			return piece.p2
		# must be last piece
		assert ((piece.p3 != a) and (piece.p3 != b))
		return piece.p3

	# out of two faces sharing an edge, return the one that is not <onePiece>
	def _getOtherFace(self, twinFaces, onePiece):
		if twinFaces[0].equals(onePiece):
			return twinFaces[1]
		assert (twinFaces[1].equals(onePiece))
		return twinFaces[0]

	# given two faces (triangles) sharing an edge <edge> return the face
	# whose 3rd point not on edge is CCW with respect to points that are on <edge>
	def _getCounterClockwiseFace(self, twinFaces, edge):
		ccwFace = twinFaces[0]
		if tools.ccw(self.points[edge[0]], self.points[edge[1]], 
					 self.points[self._getOtherPoint(ccwFace, edge[0], edge[1])]):
			return ccwFace
		return twinFaces[1]

	# draw current pieces
	def drawMe(self):
		lines = []

		for i in xrange(0, self.N):
			for j in xrange(i, self.N):
				if (self.edgeToFace[i][j][0] is not None):
					lines.append([self.points[i], self.points[j]])
				elif (self.edgeToFace[i][j][1] is not None):
					lines.append([self.points[i], self.points[j]])

		# always draw outer triangle
		lines.append([self.points[0], self.points[1]])
		lines.append([self.points[1], self.points[2]])
		lines.append([self.points[2], self.points[0]])

		lc = mc.LineCollection(lines, linewidths=1)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		ax.autoscale()
		ax.margins(0.1)
		plt.show()

	# draw current pieces with a point (the query)
	def drawMeWithPoint(self, q):
		lines = []

		for i in xrange(0, self.N):
			for j in xrange(i, self.N):
				if (self.edgeToFace[i][j][0] is not None):
					lines.append([self.points[i], self.points[j]])
				elif (self.edgeToFace[i][j][1] is not None):
					lines.append([self.points[i], self.points[j]])

		# always draw outer triangle
		lines.append([self.points[0], self.points[1]])
		lines.append([self.points[1], self.points[2]])
		lines.append([self.points[2], self.points[0]])

		lc = mc.LineCollection(lines, linewidths=1)
		fig, ax = plt.subplots()
		ax.add_collection(lc)
		ax.plot([q[0]], [q[1]], marker='o', color='k', markersize = 3)
		ax.autoscale()
		ax.margins(0.1)
		plt.show()

		return True