// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm

#include "kmeans.h"

using namespace std;

Point::Point(int id_point, vector<Real>& values, std::string name)
{
    this->id_point = id_point;
    total_values = values.size();

    for (int i = 0; i < total_values; i++)
        this->values.push_back(values[i]);

    this->name = name;
    id_cluster = -1;
}

Cluster::Cluster(int id_cluster, Point point)
{
    this->id_cluster = id_cluster;

    int total_values = point.getTotalValues();

    for (int i = 0; i < total_values; i++)
        central_values.push_back(SPN_Rand());

    points.push_back(point);
}

bool Cluster::removePoint(int id_point)
{
    int total_points = points.size();

    for (int i = 0; i < total_points; i++)
    {
        if (points[i].getID() == id_point)
        {
            points.erase(points.begin() + i);
            return true;
        }
    }
    return false;
}

int KMeans::getIDNearestCenter(const Point& point)
{
    Real sum = 0.0, min_dist;
    int id_cluster_center = 0;

    for (int i = 0; i < total_values; i++)
    {
        sum += pow(clusters[0].getCentralValue(i) - point.getValue(i), 2.0);
    }

    min_dist = sqrt(sum);

    for (int i = 1; i < K; i++)
    {
        Real dist;
        sum = 0.0;

        for (int j = 0; j < total_values; j++)
        {
            sum += pow(clusters[i].getCentralValue(j) - point.getValue(j), 2.0);
        }

        dist = sqrt(sum);

        if (dist < min_dist)
        {
            min_dist = dist;
            id_cluster_center = i;
        }
    }

    return id_cluster_center;
}

vector<Point> MatToPointVector(const Matrix_rowMajor<int>& data)
{
    vector<Point> p(data.nRows());
    for (int i = 0; i < data.nRows(); i++)
    {
        vector<Real> RealRow(data.nCols());
        for (int j=0; j<data.nCols(); j++)
            RealRow[j] = (Real)data.GetVal(i,j);
        p[i] = Point(i, RealRow, "");
    }
    return p;
}

vector<Point> MatToPointVector(const DataMatrix& dataMat)
{
    const ByteMatrix& data = dataMat.ByRow();
    vector<Point> p(data.nRows());
    for (int i = 0; i < data.nRows(); i++)
    {
        vector<Real> RealRow(data.nCols());
        for (int j=0; j<data.nCols(); j++)
            RealRow[j] = (Real)data.GetVal(i,j);
        p[i] = Point(i, RealRow, "");
    }
    return p;
}

std::vector<int> KMeans::run(const DataMatrix& data, int K, int max_iterations)
{
    vector<Point> points = MatToPointVector(data);
    this->K = K;
    total_points = data.nData();
    total_values = data.nVars();
    this->max_iterations = max_iterations;
    return run(points);
}

std::vector<int> KMeans::run(vector<Point> & points)
{
    if (K > total_points)
        printf("ERROR: K > total_points");

    cout << "KMeans run... ";
    vector<int> prohibited_indexes;

    // choose K distinct values for the centers of the clusters
    for (int i = 0; i < K; i++)
    {
        while (true)
        {
            int index_point = rand() % total_points;

            if (find(prohibited_indexes.begin(), prohibited_indexes.end(),
                     index_point) == prohibited_indexes.end())
            {
                prohibited_indexes.push_back(index_point);
                points[index_point].setCluster(i);
                Cluster cluster(i, points[index_point]);
                clusters.push_back(cluster);
                break;
            }
        }
    }

    int iter = 1;

    while (true)
    {
        bool done = true;

        // associates each point to the nearest center
        for (int i = 0; i < total_points; i++)
        {
            int id_old_cluster = points[i].getCluster();
            int id_nearest_center = getIDNearestCenter(points[i]);

            if (id_old_cluster != id_nearest_center)
            {
                if (id_old_cluster != -1)
                    clusters[id_old_cluster].removePoint(points[i].getID());

                points[i].setCluster(id_nearest_center);
                clusters[id_nearest_center].addPoint(points[i]);
                done = false;
            }
        }

        // recalculating the center of each cluster
        for (int i = 0; i < K; i++)
        {
            for (int j = 0; j < total_values; j++)
            {
                int total_points_cluster = clusters[i].getTotalPoints();
                Real sum = 0.0;

                if (total_points_cluster > 0)
                {
                    for (int p = 0; p < total_points_cluster; p++)
                        sum += clusters[i].getPoint(p).getValue(j);
                    clusters[i].setCentralValue(j, sum / total_points_cluster);
                }
            }
        }

        if (done == true || iter >= max_iterations)
        {
            cout << "Break in iteration " << iter << ".\n";
            break;
        }

        iter++;
    }

    // shows elements of clusters
    cout << "Found "<<K<<" clusters of size: ";
    for (int i = 0; i < K; i++)
    {
        int total_points_cluster = clusters[i].getTotalPoints();
        cout << total_points_cluster<< ", ";
        assert(total_points_cluster>0 && "todo");
//        cout << "Cluster " << clusters[i].getID() + 1 << " has "<< total_points_cluster<< " points"<<endl;
//        cout << "Cluster " << clusters[i].getID() + 1 << endl;
//        for (int j = 0; j < total_points_cluster; j++)
//        {
//            cout << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
//            for (int p = 0; p < total_values; p++)
//                cout << clusters[i].getPoint(j).getValue(p) << " ";
//
//            string point_name = clusters[i].getPoint(j).getName();
//
//            if (point_name != "")
//                cout << "- " << point_name;
//
//            cout << endl;
//        }


//        cout << "Cluster values: ";
//
//        for (int j = 0; j < total_values; j++)
//            cout << clusters[i].getCentralValue(j) << " ";

//        cout << "\n\n";
    }
    cout << "\n";
    //create output labels
    std::vector<int> labels(total_points);
    for (int i=0; i<total_points; i++)
        labels[i] =  points[i].getCluster();
    
    return labels;
}

