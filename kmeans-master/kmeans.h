// FOR THIS FILE ONLY THE FOLLOWING LICENSE HOLDS:
/*The MIT License (MIT)

Copyright (c) 2015 Marcos Castro de Souza

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#include "../Common.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>


class Point
{
private:
	int id_point, id_cluster;
	std::vector<Real> values;
	int total_values;
	std::string name;

public:
	Point(int id_point, std::vector<Real>& values, std::string name = "");

        Point(){}
        
	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	Real getValue(int index) const
	{
		return values[index];
	}

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(Real value)
	{
		values.push_back(value);
	}

	std::string getName()
	{
		return name;
	}
};

class Cluster
{
private:
	int id_cluster;
	std::vector<Real> central_values;
	std::vector<Point> points;

public:
	Cluster(int id_cluster, Point point);
        
	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point);

	Real getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, Real value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations;
	std::vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(const Point& point);

	std::vector<int> run(std::vector<Point> & points);
        
public:
        std::vector<int> run(const DataMatrix& data, int K, int max_iterations);
};