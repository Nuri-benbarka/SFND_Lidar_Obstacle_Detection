/* \author Aaron Brown */
// Create simple 3d highway enviroment using PCL
// for exploring self-driving car sensors

#include "sensors/lidar.h"
#include "render/render.h"
#include "processPointClouds.h"
// using templates for processPointClouds so also include .cpp to help linker
#include "processPointClouds.cpp"

std::vector<Car> initHighway(bool renderScene, pcl::visualization::PCLVisualizer::Ptr& viewer)
{

    Car egoCar( Vect3(0,0,0), Vect3(4,2,2), Color(0,1,0), "egoCar");
    Car car1( Vect3(15,0,0), Vect3(4,2,2), Color(0,0,1), "car1");
    Car car2( Vect3(8,-4,0), Vect3(4,2,2), Color(0,0,1), "car2");	
    Car car3( Vect3(-12,4,0), Vect3(4,2,2), Color(0,0,1), "car3");
  
    std::vector<Car> cars;
    cars.push_back(egoCar);
    cars.push_back(car1);
    cars.push_back(car2);
    cars.push_back(car3);

    if(renderScene)
    {
        renderHighway(viewer);
        egoCar.render(viewer);
        car1.render(viewer);
        car2.render(viewer);
        car3.render(viewer);
    }

    return cars;
}


void simpleHighway(pcl::visualization::PCLVisualizer::Ptr& viewer)
{
    // ----------------------------------------------------
    // -----Open 3D viewer and display simple highway -----
    // ----------------------------------------------------
    
    // RENDER OPTIONS
    bool renderScene = false;
    std::vector<Car> cars = initHighway(renderScene, viewer);
    
    // TODO:: Create lidar sensor
    Lidar* lidar = new Lidar(cars, 0);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = lidar->scan();
    //renderRays(viewer, lidar->position, cloud);
    //renderPointCloud(viewer,cloud,"point cloud");


    // TODO:: Create point processor
    ProcessPointClouds<pcl::PointXYZ> processor ;
    std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> Cloud = processor.SegmentPlane(cloud,100,0.2);
    //renderPointCloud(viewer,Cloud.first,"point cloud 1",Color(1,0,0));
    renderPointCloud(viewer,Cloud.second,"point cloud 2",Color(1,1,1));
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> cloudClusters = processor.Clustering(Cloud.first,1.0,3,30);

    int clusterId = 0;
    std::vector<Color> colors = {Color(1,0,0),Color(0,1,0),Color(0,0,1)};

    for(auto cluster : cloudClusters){
        std::cout  << "cluster size: " ;
        processor.numPoints(cluster);
        Box box = processor.BoundingBox(cluster);
        renderPointCloud(viewer,cluster,"Obstacle"+std::to_string(clusterId),colors[clusterId]);
        renderBox(viewer,box,clusterId);
        ++clusterId;
    }
    
  
}

void cityBlock(pcl::visualization::PCLVisualizer::Ptr& viewer,ProcessPointClouds<pcl::PointXYZI>* pointProcessorI, const pcl::PointCloud<pcl::PointXYZI>::Ptr& inputCloud)
{
  



    // Load pcd file
    //pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_raw = pointProcessorI->loadPcd("../src/sensors/data/pcd/data_1/0000000000.pcd");
    //filter cloud
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_filtered =
        pointProcessorI->FilterCloud(inputCloud,0.3,Eigen::Vector4f(-20,-8,-4,1),Eigen::Vector4f(50,5.5,8,1));
    
    // segment cloud
    std::pair<pcl::PointCloud<pcl::PointXYZI>::Ptr, pcl::PointCloud<pcl::PointXYZI>::Ptr> cloud_seg = 
        pointProcessorI->SegmentPlane(cloud_filtered,100,0.2);

    //render point cloud    
    //renderPointCloud(viewer,cloud_filtered,"inputcloud");
    renderPointCloud(viewer,cloud_seg.second,"point cloud 2",Color(1,1,1));

    //create clusters
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> cloudClusters = 
        pointProcessorI->Clustering(cloud_seg.first,0.5,5,300);

    int clusterId = 0;
    std::vector<Color> colors = {Color(1,0,0),Color(0,1,0),Color(0,0,1)};

    for(auto cluster : cloudClusters){
        std::cout  << "cluster size: " ;
        pointProcessorI->numPoints(cluster);
        Box box = pointProcessorI->BoundingBox(cluster);
        renderPointCloud(viewer,cluster,"Obstacle"+std::to_string(clusterId),colors[clusterId]);
        renderBox(viewer,box,clusterId);
        ++clusterId;
    }
    
  
}

//setAngle: SWITCH CAMERA ANGLE {XY, TopDown, Side, FPS}
void initCamera(CameraAngle setAngle, pcl::visualization::PCLVisualizer::Ptr& viewer)
{

    viewer->setBackgroundColor (0, 0, 0);
    
    // set camera position and angle
    viewer->initCameraParameters();
    // distance away in meters
    int distance = 16;
    
    switch(setAngle)
    {
        case XY : viewer->setCameraPosition(-distance, -distance, distance, 1, 1, 0); break;
        case TopDown : viewer->setCameraPosition(0, 0, distance, 1, 0, 1); break;
        case Side : viewer->setCameraPosition(0, -distance, 0, 0, 0, 1); break;
        case FPS : viewer->setCameraPosition(-10, 0, 0, 0, 0, 1);
    }

    if(setAngle!=FPS)
        viewer->addCoordinateSystem (1.0);
}


int main (int argc, char** argv)
{
    std::cout << "starting enviroment" << std::endl;

    pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    CameraAngle setAngle = XY;
    initCamera(setAngle, viewer);
    //simpleHighway(viewer);
    ProcessPointClouds<pcl::PointXYZI>* pointProcessorI = new ProcessPointClouds<pcl::PointXYZI>();
    std::vector<boost::filesystem::path> stream = pointProcessorI->streamPcd("../src/sensors/data/pcd/data_1");
    auto streamIterator = stream.begin();
    pcl::PointCloud<pcl::PointXYZI>::Ptr inputCloudI;
    //cityBlock(viewer);

    while (!viewer->wasStopped ())
    {
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

        inputCloudI = pointProcessorI->loadPcd((*streamIterator).string());
        cityBlock(viewer, pointProcessorI, inputCloudI);

        streamIterator++;
        if(streamIterator == stream.end())
            streamIterator = stream.begin();

        viewer->spinOnce ();
    } 
}