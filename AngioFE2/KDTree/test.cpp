#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "kdtree.h"
#include "catch.hpp"
#include <random>
#include <vector>
#include <chrono>

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Two dimensional test", "[KDTree]" ) {
    std::vector<double> unit(2, 1.0);
    
    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    
    kdtree.insert({-8.0, 8.0});
    kdtree.insert({-5,7});
    kdtree.insert({-2,5});
    kdtree.insert({-7,5});
    kdtree.insert({-8,2}); //current problem
    kdtree.insert({-5,2});
    kdtree.insert({-7,0});
    
    kdtree.insert({-9,-2});
    kdtree.insert({-2,-2});
    kdtree.insert({-5,-3});
    kdtree.insert({-8,-5});
    kdtree.insert({-5,-6});
    kdtree.insert({-7,-8});
    kdtree.insert({-2,-9});
    
    std::vector<double> results = kdtree.nearest({-7.0,6.0});
    //REQUIRE(results.size() == 1);
    REQUIRE(results[0] == -7.0);
    REQUIRE(results[1] == 5.0);
    
    results = kdtree.nearest({-4,1});
    //REQUIRE(results.size() == 1);
    REQUIRE(results[0] == -5.0);
    REQUIRE(results[1] == 2.0);
    
    results = kdtree.nearest({0,0});
    //REQUIRE(results.size() == 1);
    REQUIRE(results[0] == -2.0);
    REQUIRE(results[1] == -2.0);
    
    results = kdtree.nearest({-9,1});
    //REQUIRE(results.size() == 1);
    REQUIRE(results[0] == -8.0);
    REQUIRE(results[1] == 2.0);
}

TEST_CASE( "Three dimensional test", "[KDTree]" ) {
    std::vector<double> unit(3, 1.0);
    
    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    std::vector<std::vector<double>> positions;
    std::default_random_engine reng;
    std::uniform_real_distribution<double> scale(-2000, 2000);
    for(int i =0; i < 300;i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        positions.push_back(temp);
        kdtree.insert(temp);
    }
    kdtree.verifyTree();
    
    int failures = 0;
    
    for(int i =0; i < 1000; i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        std::vector<double> results = kdtree.nearest(temp);
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0] == positions[k][0]) &&
                (results[1] == positions[k][1]) &&
                (results[2] == positions[k][2])
            ){
                    REQUIRE(true);
            }
            else
            {
                
                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                        " dn:" << ndim_distance(temp, results) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results) << 
                        " temp:" << vector_printer(temp) << std::endl;
                        std::cout << std::endl;
                            failures++;
                    kdtree.nearest(temp);
                }
                REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results));
            }
            
            
        }
    }
    
    std::cout << "failures:" << failures << std::endl;
}

TEST_CASE( "Two dimensional test large", "[KDTree]" ) {
    const int dim =2;
    std::vector<double> unit(dim, 1.0);

    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    std::vector<std::vector<double>> positions;
    std::default_random_engine reng;
    std::uniform_real_distribution<double> scale(-2000, 2000);
    for(int i =0; i < 300;i++)
    {
        std::vector<double> temp;
        for(int k =0; k < dim;k++)
        {
            temp.push_back(scale(reng));
        }
        positions.push_back(temp);
        kdtree.insert(temp);
    }
    kdtree.verifyTree();
    //kdtree.PrintTree(vector_printer);
    int failures = 0;

    for(int i =0; i < 1000; i++)
    {
        std::vector<double> temp;
        for(int k =0; k < dim;k++)
        {
            temp.push_back(scale(reng));
        }
        std::vector<double> results = kdtree.nearest(temp);
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0] == positions[k][0]) &&
                (results[1] == positions[k][1])
              ){
                  REQUIRE(true);
                }
            else
            {

                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                    " dn:" << ndim_distance(temp, results) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results) << 
                    " temp:" << vector_printer(temp) << std::endl;
                    std::cout << std::endl;
                    failures++;
                    kdtree.nearest(temp);
                }
                REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results));
            }


        }
    }

    std::cout << "failures:" << failures << std::endl;
  }
  
TEST_CASE( "Three dimensional test timed", "[KDTree]" ) {
    std::vector<double> unit(3, 1.0);

    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    std::vector<std::vector<double>> positions;
    std::default_random_engine reng;
    std::uniform_real_distribution<double> scale(-2000, 2000);
    for(int i =0; i < 1000;i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        positions.push_back(temp);
        kdtree.insert(temp);
    }
    kdtree.verifyTree();

    int failures = 0;
    
    const int num_tests = 1000;
    std::vector<std::vector<double> > testData;
    for(int i =0; i < num_tests;i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        testData.push_back(temp);
    }
    
    std::vector<std::vector<double> > resultData;
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    for(int i =0; i < num_tests;i++)
    {
        std::vector<double> results = kdtree.nearest(testData[i]);
        resultData.push_back(results);
    }
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::cout <<  "querys on the kd tree took: " <<  time_span.count() <<  std::endl;
    
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    int best =0;
    for(int i =0; i < num_tests; i++)
    {
        std::vector<double> temp = testData[i]; 
        std::vector<double> results = resultData[i];
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0] == positions[k][0]) &&
                (results[1] == positions[k][1]) &&
                (results[2] == positions[k][2])
              ){
                  //REQUIRE(true);
                }
            else
            {

                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                    " dn:" << ndim_distance(temp, results) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results) << 
                    " temp:" << vector_printer(temp) << std::endl;
                    std::cout << std::endl;
                    failures++;
                    kdtree.nearest(temp);
                    best = i;
                }
                //REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results));
            }


        }
    }
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();

    std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t4-t3);
    std::cout <<  "brute force takes: " <<  time_span2.count() <<  std::endl;
    
    for(int i =0; i < num_tests; i++)
    {
        std::vector<double> temp = testData[i]; 
        std::vector<double> results = resultData[i];
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0] == positions[k][0]) &&
                (results[1] == positions[k][1]) &&
                (results[2] == positions[k][2])
              ){
                  REQUIRE(true);
                }
            else
            {

                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                    " dn:" << ndim_distance(temp, results) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results) << 
                    " temp:" << vector_printer(temp) << std::endl;
                    std::cout << std::endl;
                    failures++;
                    kdtree.nearest(temp);
                }
                REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results));
            }


        }
    }
    
    std::cout << "failures:" << failures << std::endl;
  }

  
  
TEST_CASE( "Three dimensional test multi insert", "[KDTree]" ) {
    std::vector<double> unit(3, 1.0);

    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    std::vector<std::vector<double>> positions;
    std::default_random_engine reng;
    std::uniform_real_distribution<double> scale(-2000, 2000);
    for(int i =0; i < 30;i++)
    {
        
        std::list<std::vector<double>> mis;
        for (int j = 0; j < 8;j++)
        {
            std::vector<double> temp;
            for(int k =0; k < 3;k++)
            {
                temp.push_back(scale(reng));
            }
            positions.push_back(temp);
            //kdtree.insert(temp);
            mis.push_front(temp);
        }
        
        kdtree.insert(mis);
    }
    kdtree.verifyTree();

    int failures = 0;

    for(int i =0; i < 1000; i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        std::vector<double> results = kdtree.nearest(temp);
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0] == positions[k][0]) &&
                (results[1] == positions[k][1]) &&
                (results[2] == positions[k][2])
              ){
                  REQUIRE(true);
                }
            else
            {

                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                    " dn:" << ndim_distance(temp, results) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results) << 
                    " temp:" << vector_printer(temp) << std::endl;
                    std::cout << std::endl;
                    failures++;
                    kdtree.nearest(temp);
                }
                REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results));
            }


        }
    }

    std::cout << "failures:" << failures << std::endl;
  }
  
  /*
TEST_CASE( "Three dimensional test multi insert within", "[KDTree]" ) {
    std::vector<double> unit(3, 1.0);

    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    std::vector<std::vector<double>> positions;
    std::default_random_engine reng;
    std::uniform_real_distribution<double> scale(-2000, 2000);
    for(int i =0; i < 300;i++)
    {

        std::list<std::vector<double>> mis;
        for (int j = 0; j < 8;j++)
        {
            std::vector<double> temp;
            for(int k =0; k < 3;k++)
            {
                temp.push_back(scale(reng));
            }
            positions.push_back(temp);
            //kdtree.insert(temp);
            mis.push_front(temp);
        }

        kdtree.insert(mis);
    }
    kdtree.verifyTree();

    int failures = 0;
    bool ever_in_bounds = false;
    for(int i =0; i < 1000; i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        double within_d = 100.0*100.0;
        std::vector<std::vector<double>> results = kdtree.within(temp,  within_d);
        for(int k=0; k < positions.size(); k++)
        {
            bool wid = ndim_distance(temp,  positions[k]) <= within_d; 
            bool found = false;
            for(int q =0; q < results.size();q++)
            {
                if(results[q][0] == positions[k][0] &&
                    results[q][1] == positions[k][1] &&
                    results[q][2] == positions[k][2])
                    {
                        found = true;
                        ever_in_bounds = true;
                        break;
                    }
            }
            REQUIRE(found ==  wid);
        }
    }
    REQUIRE(ever_in_bounds ==  true);

    std::cout << "failures:" << failures << std::endl;
  }
 
  
TEST_CASE( "Three dimensional test multi insert rebuild", "[KDTree]" ) {
    std::vector<double> unit(3, 1.0);
    
    
    auto kdtree = KDTree<std::vector<double>, std::vector<double>>(null_accessor, ndim_distance, ndim_distance_to_plane, unit);
    std::vector<std::vector<double>> positions;
    std::default_random_engine reng;
    std::uniform_real_distribution<double> scale(-2000, 2000);
    for(int i =0; i < 30;i++)
    {

        std::list<std::vector<double>> mis;
        for (int j = 0; j < 8;j++)
        {
            std::vector<double> temp;
            for(int k =0; k < 3;k++)
            {
                temp.push_back(scale(reng));
            }
            positions.push_back(temp);
            //kdtree.insert(temp);
            mis.push_front(temp);
        }

        kdtree.insert(mis);
    }
    kdtree.verifyTree();

    int failures = 0;

    for(int i =0; i < 1000; i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        std::vector<std::vector<double>> results = kdtree.nearest(temp);
        REQUIRE(results.size() == 1);
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0][0] == positions[k][0]) &&
                (results[0][1] == positions[k][1]) &&
                (results[0][2] == positions[k][2])
              ){
                  REQUIRE(true);
                }
            else
            {

                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results[0]))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                    " dn:" << ndim_distance(temp, results[0]) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results[0]) << 
                    " temp:" << vector_printer(temp) << std::endl;
                    std::cout << std::endl;
                    failures++;
                    kdtree.nearest(temp);
                }
                REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results[0]));
            }


        }
    }
    kdtree.rebuild();
    kdtree.rebuild();
    kdtree.rebuild();
    kdtree.rebuild();
    
    for(int i =0; i < 1000; i++)
    {
        std::vector<double> temp;
        for(int k =0; k < 3;k++)
        {
            temp.push_back(scale(reng));
        }
        std::vector<std::vector<double>> results = kdtree.nearest(temp);
        REQUIRE(results.size() == 1);
        for(int k=0; k < positions.size(); k++)
        {
            if((results[0][0] == positions[k][0]) &&
                (results[0][1] == positions[k][1]) &&
                (results[0][2] == positions[k][2])
              ){
                  REQUIRE(true);
                }
            else
            {

                if(ndim_distance(positions[k], temp) < ndim_distance(temp, results[0]))
                {
                    std::cout << "dk:" << ndim_distance(positions[k], temp) <<
                    " dn:" << ndim_distance(temp, results[0]) << std::endl;
                    std::cout  << "pos: " << vector_printer(positions[k]) << " res:" << vector_printer( results[0]) << 
                    " temp:" << vector_printer(temp) << std::endl;
                    std::cout << std::endl;
                    failures++;
                    kdtree.nearest(temp);
                }
                REQUIRE(ndim_distance(positions[k], temp) >= ndim_distance(temp, results[0]));
            }


        }
    }

    std::cout << "failures:" << failures << std::endl;
  }
 */
