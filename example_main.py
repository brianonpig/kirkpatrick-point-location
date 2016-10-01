from random import randint
import itertools
from Kirkpatrick import Kirkpatrick

def main():

	# randomly generate <numpoints> points and perform point location
	numPoints = 50
	points = [[randint(0, 100000), randint(0, 100000)] for i in xrange(numPoints)]
	points.sort()
	points = list(points for points,_ in itertools.groupby(points))

	kp = Kirkpatrick(points)
	query = [randint(0,100000), randint(0,100000)]
	print kp.animatedLocation(query)

if __name__ == "__main__":
    main()