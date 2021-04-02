/*
        PSOVina version 2.0  

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
#ifndef PSO_H_
#define PSO_H_

#include "common.h"
#include "conf.h"
#include "model.h"
#include <vector>

class pso
{	
public:
	struct bird{
		     vec velocity,vO;  //velocity of position, velocity of orientation
		     
		     vec pbest_pos;
		     vec current_position;
 
		     qt pbest_orientation; 
		     qt current_orientation;
		     
		     fl* pbest_torsion;
		     fl* current_torsion;
		     fl* vT;    //velocity of torsion

             int order;
		     double pbest_fit;  
             double current_fit;

	};
	
	int torsionSize;
	double c1,c2;						// weight, learning coefficient(1&2)
	rng g;

	int number;							//number of birds
	vec corner1,corner2;				//corners of search box
	
	vec gbest_position;			//global best of degree of freedom vector
	qt gbest_orientation;		//global best of degree of freedom vector
	fl* gbest_torsion;			//global best of degree of freedom vector
	double gbest_fit;			//global best value

	std::vector<bird> particle;
	
	
	  double R1Max_;
	  double R1Min_;
	  double R2Max_;
	  double R2Min_;
	
	
	pso(int,double,double,const vec,const vec,rng&,conf&);
	void init(rng&,conf&);

    void computeQPSOPositions(int, fl);
    void computeQPSOOrientation(int, fl);
    void computeQPSOTorsion(int, int, fl);

    /*Update Velocity for vector*/
    void updateVelocity(rng&,int,double,double);
    void updateVelocityO(rng&,int,double,double);
    void updateVelocityT(rng&,int,sz,double,double);

    /*Compute the new pose*/
    void computePSOPositions(int);
    void computePSOOrientation(int);
    void computePSOTorsion(int,rng&,sz);

    void computeRDPSOPositions(int, fl, fl);
    void computeRDPSOOrientation(int, fl, fl);
    void computeRDPSOTorsion(int,sz, fl, fl);

    void computePositions(int, fl);
    void computeOrientation(int, fl);
    void computeTorsion(int,sz, fl);

    vec calculateAveragePosition();
    qt calculateAverageOrientation();
    fl calculateAverageTorsion(sz);

    vec calculateAveragePosition1();
    qt calculateAverageOrientation1();
    fl calculateAverageTorsion1(sz);

    /*Update personal best value*/
	void updatePersonalBest(int,double);
    void updateCurrentBest(int,double);
    
    /*Update global best value*/
	void updateGlobalBest(int);
	void updateGlobalBest_1(fl);
    
    /*Get personal best value*/
	double getPersonalBest(int);	//return the personal best value
    double getCurrentBest(int);    //return the personal best value

    /*Set personal best vector*/
    void updateBestPosition(int,vec);
    void updateBestOrientation(int,qt);
    void updateBestTorsion(int,fl,sz);
    void updateCurrentPosition(int,vec);
    void updateCurrentOrientation(int,qt);
    void updateCurrentTorsion(int,fl,sz);

    /*Get current vector*/
	vec getCurrentPosition(int);
    qt getCurrentOrientation(int);
    fl getCurrentTorsion(int,sz);
    vec getPBestPosition(int);
    qt getPBestOrientation(int);
    fl getPBestTorsion(int,sz);
};


#endif /*PSO_H_*/
