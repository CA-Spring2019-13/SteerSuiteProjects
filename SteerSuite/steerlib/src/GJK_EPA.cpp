#include "obstacles/GJK_EPA.h"

// Vector A = Support(shape A, init) - Support (Shape B, - init)
// Simplex S = {A}
// Vector D = -A
/*
	A = Support(Shape A, D) - Support(Shape B, -D)
	if dot prod < 0 --> false


*/

SteerLib::GJK_EPA::GJK_EPA()
{
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	//Initialize a simplex
	std::vector<Util::Vector> simplex; 

	//Minkowski Center Calculation
	Util::Vector centerA = locateCenter(_shapeA);
	Util::Vector centerB = locateCenter(_shapeB);

	//If the centers are the same, there is definitely a collision
	if (centerA == centerB) {

		EPA(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB, simplex);
		return true; // A collision has occured!
	}

	//If the above condition didnt occur, begin adding points to our simplex vector
	Util::Vector s1 = support(_shapeA, _shapeB, (centerA - centerB));
	Util::Vector s2 = support(_shapeA, _shapeB, (-1 * (centerA - centerB)));
	simplex.push_back(s1);
	simplex.push_back(s2);
	Util::Vector ab;
	Util::Vector ao;


	//Continue to add points to our simplex vector, return if we find a collision
	while (true) {

		int size = simplex.size();
		Util::Vector o;
		o.x = 0;
		o.y = 0;
		o.z = 0;
		//Minkowski diff for A & B, and A & Origin
		ab = simplex[size - 1] - simplex[size - 2];
		ao = o - simplex[size - 1];

		//Cross product to determine the new direction to traverse (find next max outlier)
		//Use Lagrange formula and/or Jacobi identity
		Util::Vector D = ao * (ab*ab);
		D = D - (ab*(ab*ao));

		//check if origin

		//Add to our simplex after sending to support function
		Util::Vector s = support(_shapeA, _shapeB, D);
		simplex.push_back(s);

		//Check if the simplex now contains the origin
		if (simplex[2] * D > 0) {

			if (hasOrigin(simplex)) {
				EPA(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB, simplex);
				return true; //Contains origin, collision found!
			}

		}
		else {
			//Not in given area
			return false;
		}




		






	}






	
	
	
	
	return false; // There is no collision
}


Util::Vector locateCenter(const std::vector<Util::Vector>& myShape) {

	//Create a vector for the center
	Util::Vector cent(0, 0, 0);

	//If the size is 0, return the default center (can't divide by zero for avg calc)
	int size = myShape.size();
	if (size == 0) {
		return cent;
	}

	int xTotal = 0;
	int yTotal = 0;
	int zTotal = 0;
	int xAvg, yAvg, zAvg;
	int i;

	//Find average in x direction
	for (i = 0; i < size; i++) {

		xTotal += myShape[i].x;
	}

	xAvg = (xTotal) / size;

	//Find average in y direction

	for (i = 0; i < size; i++) {

		yTotal += myShape[i].y;
	}

	yAvg = (yTotal) / size;

	//Find average in z direction

	for (i = 0; i < size; i++) {

		zTotal += myShape[i].z;
	}

	zAvg = (zTotal) / size;

	cent.x = xAvg;
	cent.y = yAvg;
	cent.z = zAvg;

	return cent;
}

Util::Vector maxDistanceInD(Util::Vector D, const std::vector<Util::Vector>& myShape)
{

	//Initialize a vector to begin tracking distances
	Util::Vector outlier(0, 0, 0);
	float max = 0;
	int size = myShape.size();

	for (int i = 0; i < size; i++) {

		if (((myShape[i])* D) > max) {
			outlier = myShape[i];
			//Search for the largest dot product in the D direction
			//Will be used for subtraction in support
			max = myShape[i] * D;



		}
	}

	return outlier;
}

Util::Vector support(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, Util::Vector D)
{
	/*Support function:
	1. Furthest point of A in D direction
	2. Furthest point of B in -D direction
	3. Subtract #1 from #2, return that result */

	Util::Vector negD = -1 * D;
	return (maxDistanceInD(D, _shapeA) - maxDistanceInD(negD, _shapeB));

}

void EPA(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector> &simplex) {

	float currentDistance;

	Util::Vector a, b, ab, abao, V, ao;
	int size = simplex.size();
	float magAB;
	float magAO;
	int loc = 0;
	float shortestDistance = FLT_MAX;
	float closestPointDistance;


	//Calculate the penetration
	while (true) {

		for (int i = 0; i < size; i++) {

			a = simplex[i];

			if (i == (size - 1)) {
				b = simplex[0];
			}
			else {
				b = simplex[i + 1];
			}

			//Difference, then calc magnitudes:
			ab = b - a;
			ao = -1 * a;
			magAB = ab.norm();
			magAO = ao.norm();
			//abao = fabs(ab * ao);

			//Use absolute value for distances fabs()
			//Magnitude of the vector vect.norm()
			if (fabs(fabs(ab*ao) - ab.norm()*ao.norm()) >= FLT_MIN) {

				V = ao * (ab*ab) - ab * (ab*ao);
				V = Util::normalize(V);

			}
			else {
				Util::Vector normVect(-1 * ab.z, 0, ab.x);
				V = Util::normalize(normVect);
			}


			//Calculate the distance
			currentDistance = fabs(ao*V);
			//distance = fabs(distance);

			//Normalize again
			V = Util::normalize(V);

			//Now check if the current distance is the smallest distance
			//Initialize to max float, then reduce it whenever a smaller dist is found


			if (currentDistance < shortestDistance) {

				V = -1 * V;

				if (i == size - 1) {

					loc = 0;
				}
				else {

					loc = i + 1;
				}

				shortestDistance = currentDistance;








				//check that float is greater than 0
				//instead of 0 use FLT_MIN





			}




		}

		Util::Vector s = support(_shapeA, _shapeB, V);
		closestPointDistance = s * V - shortestDistance;
		closestPointDistance = fabs(closestPointDistance);

		if (closestPointDistance > FLT_MIN) {

			simplex.insert(simplex.begin() + loc, s);
		}
		else {

			return_penetration_depth = shortestDistance;
			return_penetration_vector = V;
			return;
		}






	}
}



	bool hasOrigin(std::vector<Util::Vector> &simplex) {

		Util::Vector a, b, c, aNeg, ab, ac;

		//Need to determine if the origin lies in the area formed by the simplex
		//Cross product = area of the parallelogram between two vectors
		//If 0 or negative between perpendicular then we can say origin is not in the area
		//Otherwise return true.

		//Third vector
		a = simplex[2];
		aNeg = -1 * a;
		//Second vector
		b = simplex[1];
		//First vector
		c = simplex[0];

		//Cross product of given vectors
		ab = (((c - a)*((b - a)*(b - a))) - ((b - a)*((b - a)*(c - a))));
		ac = (((b - a)*((c - a)*(c - a))) - ((c - a)*((c - a)*(b - a))));

		if ((ab * (-1 * a)) <= 0) {
			//Get rid of the point because it is no longer of use
			simplex.erase(simplex.begin());
			return false;
		}

		if ((ac * (-1 * a)) <= 0) {
			simplex.erase(simplex.begin() + 1);
			return false;
		}


		return true; //If the above conditions fail to hold, then the given area has the origin








	}






	
















