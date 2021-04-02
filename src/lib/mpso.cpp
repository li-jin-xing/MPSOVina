/*
        PSOVina version 2.0 

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
#include "mpso.h"
#include "random.h"

	pso::pso(int num_birds,double c1,double c2,const vec corner1,const vec corner2, rng& g,conf& c)
	{
		  sz torsionSize = c.ligands[0].torsions.size();
		  this->c1 = c1;
		  this->c2 = c2;
		  this->number = num_birds;
		  this->g = g;
		  this->corner1[0] = corner1[0];			//minmum
		  this->corner1[1] = corner1[1];
		  this->corner1[2] = corner1[2];
		  this->corner2[0] = corner2[0];			//maximum
		  this->corner2[1] = corner2[1];
		  this->corner2[2] = corner2[2];
		  this->torsionSize = (int)torsionSize;
		  this->R1Max_ = 1;
		  this->R1Min_ = 0;
		  this->R2Max_ = 1;
		  this->R2Min_ = 0;
		  this->gbest_torsion = new fl[torsionSize];
		  init(g,c);
	}
	
	void pso::init(rng &g,conf& c)
	{

		int i;
		for(i=0;i<this->number;i++)
		{
			bird single_bird;

			single_bird.pbest_fit = 1.7976931348623158e+308;
			single_bird.current_fit = 1.7976931348623158e+308;
			single_bird.order=0;
			//set position part
			single_bird.velocity = random_in_box(this->corner1,this->corner2,g);
			single_bird.current_position = random_in_box(this->corner1,this->corner2,g);
			
			//set orientation part
			single_bird.vO = random_inside_sphere(g);
			qt tmp_o = c.ligands[0].rigid.orientation;
			quaternion_increment(tmp_o,  random_inside_sphere(g));
			single_bird.current_orientation = tmp_o;
			
			//init. the array for the number of torsion
			single_bird.current_torsion=new fl[this->torsionSize];
			single_bird.vT=new fl[this->torsionSize];
			single_bird.pbest_torsion=new fl[this->torsionSize];
			
			for(int x=0;x<this->torsionSize;x++)						//init. all the torsion that the ligand has
			{
				single_bird.vT[x] = random_fl(-pi, pi, g);
				single_bird.current_torsion[x] = random_fl(-pi, pi, g);
			
			}
			
			particle.push_back(single_bird);
		}
		
		this->gbest_fit = 1.7976931348623158e+308;
	}

    void pso::computeQPSOPositions(int i, fl a)
    {
        vec ave_position = this->calculateAveragePosition1();
        for (int k=0; k<3; k++)
        {
            fl fi1, fi2, p, v, b, z;
            fi1 = random_fl(0,1,this->g);
            fi2 = random_fl(0,1,this->g);
            p = (fi1*this->c1*particle[i].pbest_pos[k] + fi2*this->c2*this->gbest_position[k])/(fi1*this->c1+fi2*this->c2);
            v = std::log(1.0/random_fl(0,1,this->g));
            b = a * std::fabs(ave_position[k]-particle[i].current_position[k]);

            z = random_fl(0,1,this->g);
            if (z < 0.5)
                particle[i].current_position[k] = p + b * v;
            else
                particle[i].current_position[k] = p - b * v;
            
        }
        //give a random position, if outside the search box
        if(particle[i].current_position[0] < corner1[0])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[1] < corner1[1])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[2] < corner1[2])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[0] > corner2[0])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[1] > corner2[1])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[2] > corner2[2])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
    }

    void pso::computeQPSOOrientation(int i, fl a)
    {
        qt ave_orientation = this->calculateAverageOrientation1(), p, b;
        fl z,v;
        qt fi1 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
        qt fi2 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
        p = quaternion_div(quaternion_mul(fi1*this->c1, particle[i].pbest_orientation)+quaternion_mul(fi2*this->c2, this->gbest_orientation), (fi1*this->c1+fi2*this->c2));
        b = ave_orientation - particle[i].current_orientation;
        v = std::log(1.0/random_fl(0,1,this->g));

        z = random_fl(0,1,this->g);
        if (z < 0.5)
            quaternion_increment(p, a*v*vec_abs(quaternion_to_angle(b)));
        else
            quaternion_increment(p, -a*v*vec_abs(quaternion_to_angle(b)));
    }

    void pso::computeQPSOTorsion(int i ,int which, fl a)
    {
        fl ave_torsion = this->calculateAverageTorsion1(which);
        
        fl fi1, fi2, p, v, b, z;
        fi1 = random_fl(0,1,this->g);
        fi2 = random_fl(0,1,this->g);
        p = (fi1*this->c1*particle[i].pbest_torsion[which] + fi2*this->c2*this->gbest_torsion[which])/(fi1*this->c1+fi2*this->c2);
        v = std::log(1.0/random_fl(0,1,this->g));
        b = a * std::fabs(ave_torsion-particle[i].current_torsion[which]);

        z = random_fl(0,1,this->g);
        if (z < 0.5)
            particle[i].current_torsion[which] = p + b * v;
        else
            particle[i].current_torsion[which] = p - b * v;

        particle[i].current_torsion[which] = normalized_angle(particle[i].current_torsion[which]);
        if(isnan(particle[i].current_torsion[which]))
            particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
        if(isinf(particle[i].current_torsion[which]))
            particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
    }

    void pso::updateVelocity(rng& generator,int i,double cm,double l)
    {
        particle[i].velocity[0] = particle[i].velocity[0]*l+this->c1*cm*(particle[i].pbest_pos[0]-particle[i].current_position[0])+this->c2*(1-cm)*(this->gbest_position[0]-particle[i].current_position[0]);
        particle[i].velocity[1] = particle[i].velocity[1]*l+this->c1*cm*(particle[i].pbest_pos[1]-particle[i].current_position[1])+this->c2*(1-cm)*(this->gbest_position[1]-particle[i].current_position[1]);
        particle[i].velocity[2] = particle[i].velocity[2]*l+this->c1*cm*(particle[i].pbest_pos[2]-particle[i].current_position[2])+this->c2*(1-cm)*(this->gbest_position[2]-particle[i].current_position[2]);
    }
    
    void pso::updateVelocityO(rng& generator,int i,double cm,double l)
    {
        qt p1 = particle[i].pbest_orientation-particle[i].current_orientation;
        qt p2 = this->gbest_orientation-particle[i].current_orientation;
        particle[i].vO = particle[i].vO*l+this->c1*cm*quaternion_to_angle(p1)+this->c2*(1-cm)*quaternion_to_angle(p2);
    }
    
    void pso::updateVelocityT(rng& generator,int i,sz which,double cm,double l)
    {
        particle[i].vT[which] = particle[i].vT[which]*l+this->c1*cm*(particle[i].pbest_torsion[which]-particle[i].current_torsion[which])+this->c2*(1-cm)*(this->gbest_torsion[which]-particle[i].current_torsion[which]);
    }

    void pso::computePSOPositions(int i)
    {
        particle[i].current_position[0] = particle[i].current_position[0] + particle[i].velocity[0];
        particle[i].current_position[1] = particle[i].current_position[1] + particle[i].velocity[1];
        particle[i].current_position[2] = particle[i].current_position[2] + particle[i].velocity[2];

        //give a random position, if outside the search box
        if(particle[i].current_position[0] < corner1[0])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[1] < corner1[1])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[2] < corner1[2])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        
        if(particle[i].current_position[0] > corner2[0])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[1] > corner2[1])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[2] > corner2[2])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
    }


    void pso::computeRDPSOPositions(int i, fl alpha, fl beta)
    {
        vec ave_position = this->calculateAveragePosition();
        fl fi1, fi2, p, v, rdn, V;
        for (int k=0; k<3; k++)
        {
            fi1 = random_fl(0,1,this->g);
            fi2 = random_fl(0,1,this->g);
            p = (fi1*particle[i].pbest_pos[k] + fi2*pso::gbest_position[k])/(fi1+fi2);
            v = std::log(1.0/random_fl(0,1,this->g));
            rdn = random_fl(0,1,this->g);
            V=alpha*fabs(ave_position[k]-particle[i].current_position[k])*rdn+beta*(p-particle[i].current_position[k]);

            if (V>(corner2[k]-corner1[k]))
                V=corner2[k]-corner1[k];
            if (V<(corner1[k]-corner2[k]))
                V=corner1[k]-corner2[k];
            particle[i].current_position[k] = particle[i].current_position[k] + V;
        }

        //give a random position, if outside the search box
        if(particle[i].current_position[0] < corner1[0])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[1] < corner1[1])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[2] < corner1[2])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        
        if(particle[i].current_position[0] > corner2[0])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[1] > corner2[1])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
        if(particle[i].current_position[2] > corner2[2])
            particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
    }

    void pso::computeRDPSOOrientation(int i, fl alpha, fl beta)
    {
        qt ave_orientation = this->calculateAverageOrientation(), p, p1, p2;
        fl rdn;

        qt fi1 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
        qt fi2 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
        p = quaternion_div(quaternion_mul(fi1, particle[i].pbest_orientation)+quaternion_mul(fi2, pso::gbest_orientation), (fi1+fi2));
        rdn=random_fl(0,1,this->g);

        p1=ave_orientation - particle[i].current_orientation;
        p2=p-particle[i].current_orientation;

        vec tmp_v = alpha*vec_abs(quaternion_to_angle(p1))*rdn + beta*quaternion_to_angle(p2);
        quaternion_increment(particle[i].current_orientation, tmp_v);
    }

    void pso::computeRDPSOTorsion(int i,sz which, fl alpha, fl beta)
    {
        fl ave_torsion = this->calculateAverageTorsion(which);
        fl fi1, fi2, p, rdn, v, V;
        fi1 = random_fl(0,1,this->g);
        fi2 = random_fl(0,1,this->g);
        p = (fi1*particle[i].pbest_torsion[which] + fi2*pso::gbest_torsion[which])/(fi1+fi2);
        v = std::log(1.0/random_fl(0,1,this->g));
        rdn=random_fl(0,1,this->g);
        V=alpha*std::fabs(ave_torsion-particle[i].current_torsion[which])*rdn+beta*(p-particle[i].current_torsion[which]);
        particle[i].current_torsion[which] = particle[i].current_torsion[which] + V;

        particle[i].current_torsion[which] = normalized_angle(particle[i].current_torsion[which]);
        if(isnan(particle[i].current_torsion[which]))
            particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
        if(isinf(particle[i].current_torsion[which]))
            particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
    }

    void pso::computePSOOrientation(int i)
    {
        vec tmp_v = particle[i].vO;
        quaternion_increment(particle[i].current_orientation, tmp_v);
    }

    void pso::computePSOTorsion(int i,rng& generator,sz which)
    {
        particle[i].current_torsion[which] = particle[i].current_torsion[which] + particle[i].vT[which];
        particle[i].current_torsion[which] = normalized_angle(particle[i].current_torsion[which]);

        if(isnan(particle[i].current_torsion[which]))
            particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
        if(isinf(particle[i].current_torsion[which]))
            particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
        particle[i].current_torsion[which] = normalized_angle(particle[i].current_torsion[which]);
    }
    
	void pso::computePositions(int i, fl alpha)
	{
		fl fi1, fi2, p;
		for (int k=0; k<3; k++)
		{
			fi1 = random_fl(0,1,this->g);
			fi2 = random_fl(0,1,this->g);
			p = (fi1*this->c1*particle[i].pbest_pos[k] + fi2*this->c2*this->gbest_position[k])/(fi1*this->c1+fi2*this->c2);
			particle[i].current_position[k] = particle[i].current_position[k] + alpha*(p-particle[i].current_position[k]);
		}

		//give a random position, if outside the search box
		if(particle[i].current_position[0] < corner1[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] < corner1[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] < corner1[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		
		if(particle[i].current_position[0] > corner2[0])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[1] > corner2[1])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
		if(particle[i].current_position[2] > corner2[2])
			particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
	}

	void pso::computeOrientation(int i, fl alpha)
	{
		qt fi_1 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt fi_2 = qt(random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g), random_fl(0,1,this->g));
		qt p = quaternion_div(quaternion_mul(fi_1*this->c1, particle[i].pbest_orientation)+quaternion_mul(fi_2*this->c2, this->gbest_orientation), (fi_1*this->c1+fi_2*this->c2));
		qt fi2 = p-particle[i].current_orientation;
		quaternion_increment(particle[i].current_orientation, alpha*quaternion_to_angle(fi2));
	}

	void pso::computeTorsion(int i,sz which, fl alpha)
	{
		fl fi1, fi2, p;
		fi1 = random_fl(0,1,this->g);
		fi2 = random_fl(0,1,this->g);
		p = (fi1*this->c1*particle[i].pbest_torsion[which] + fi2*this->c2*this->gbest_torsion[which])/(fi1*this->c1+fi2*this->c2);
		particle[i].current_torsion[which] = particle[i].current_torsion[which] + alpha*(p-particle[i].current_torsion[which]);

		particle[i].current_torsion[which] = normalized_angle(particle[i].current_torsion[which]);
		if(isnan(particle[i].current_torsion[which]))
			particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
		if(isinf(particle[i].current_torsion[which]))
			particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
	}

	void pso::updatePersonalBest(int i,double e)
	{
		particle[i].pbest_fit = e;
	}

	void pso::updateCurrentBest(int i,double e)
	{
		particle[i].current_fit = e;
	}
	
	void pso::updateGlobalBest(int i)
	{
		this->gbest_fit = particle[i].pbest_fit;
	}
	
	void pso::updateGlobalBest_1(fl en)
	{
		this->gbest_fit = en;
	}

	
	double pso::getPersonalBest(int i)
	{
		return particle[i].pbest_fit;
	}

	void pso::updateBestPosition(int i,vec pos)
	{
		particle[i].pbest_pos = pos;
		
	}
	
	void pso::updateBestOrientation(int i, qt orientation)
	{
		particle[i].pbest_orientation = orientation;
	}
	
	void pso::updateBestTorsion(int i, fl torsion,sz which)
	{
		particle[i].pbest_torsion[which] = torsion;
	}

	void pso::updateCurrentPosition(int i,vec pos)
	{
		particle[i].current_position = pos;
		
	}
	
	void pso::updateCurrentOrientation(int i, qt orientation)
	{
		particle[i].current_orientation = orientation;
	}
	
	void pso::updateCurrentTorsion(int i, fl torsion,sz which)
	{
		particle[i].current_torsion[which] = torsion;
	}
	
	vec pso::getCurrentPosition(int i)
	{
		return particle[i].current_position;
	}
	
	qt pso::getCurrentOrientation(int i)
	{
		return particle[i].current_orientation;
	}
	
	fl pso::getCurrentTorsion(int i,sz which)
	{
		return particle[i].current_torsion[which];
	}

	fl pso::getCurrentBest(int i)
	{
		return particle[i].current_fit;
	}

	vec pso::getPBestPosition(int i)
	{
		return particle[i].pbest_pos;
	}
	
	qt pso::getPBestOrientation(int i)
	{
		return particle[i].pbest_orientation;
	}
	
	fl pso::getPBestTorsion(int i,sz which)
	{
		return particle[i].pbest_torsion[which];
	}

	vec pso::calculateAveragePosition()
	{
		vec tmp = vec(0, 0, 0);
		for(int i=0;i<this->number;i++)
		{
			tmp[0] += particle[i].pbest_pos[0];
			tmp[1] += particle[i].pbest_pos[1];
			tmp[2] += particle[i].pbest_pos[2];
		}
		tmp[0] = tmp[0] / float(this->number);
		tmp[1] = tmp[1] / float(this->number);
		tmp[2] = tmp[2] / float(this->number);
		return tmp;
	}

	qt pso::calculateAverageOrientation()
	{
		float tmp[4] = { 0, 0, 0, 0 };
		float len = this->number;
		for (int i = 0; i < len; i++)
		{
			tmp[0] += particle[i].pbest_orientation.R_component_1();
			tmp[1] += particle[i].pbest_orientation.R_component_2();
			tmp[2] += particle[i].pbest_orientation.R_component_3();
			tmp[3] += particle[i].pbest_orientation.R_component_4();
		}
		return qt(1.0/len*tmp[0], 1.0/len*tmp[1], 1.0/len*tmp[2], 1.0/len*tmp[3]);
	}

	fl pso::calculateAverageTorsion(sz which)
	{
		fl tmp = 0;
		for(int i=0; i<this->number; i++)
			tmp += particle[i].pbest_torsion[which];
		tmp = tmp / float(this->number);
		return tmp;
	}


    vec pso::calculateAveragePosition1()
    {
        vec tmp = vec(0, 0, 0);
        for(int i=0;i<this->number;i++)
        {
            tmp[0] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_pos[0];
            tmp[1] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_pos[1];
            tmp[2] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_pos[2];
        }
        tmp[0] = tmp[0] / float(this->number);
        tmp[1] = tmp[1] / float(this->number);
        tmp[2] = tmp[2] / float(this->number);
        return tmp;
    }

    qt pso::calculateAverageOrientation1()
    {
        float tmp[4] = { 0, 0, 0, 0 };
        float len = this->number;
        for (int i = 0; i < len; i++)
        {
            tmp[0] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_orientation.R_component_1();
            tmp[1] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_orientation.R_component_2();
            tmp[2] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_orientation.R_component_3();
            tmp[3] += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_orientation.R_component_4();
        }
        return qt(1.0/len*tmp[0], 1.0/len*tmp[1], 1.0/len*tmp[2], 1.0/len*tmp[3]);
    }

    fl pso::calculateAverageTorsion1(sz which)
    {
        fl tmp = 0;
        for(int i=0; i<this->number; i++)
            tmp += particle[i].pbest_fit/pso::gbest_fit*particle[i].pbest_torsion[which];
        tmp = tmp / float(this->number);
        return tmp;
    }