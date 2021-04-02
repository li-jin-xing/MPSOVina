/*
       PSOVina version 2.0  

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Giotto H. K. TAI  <giottotai@yahoo.com.hk>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "mpso_mutate.h"
#include <math.h>
#include <time.h>
#define PI 3.14159265

sz pso_count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

void mpso_mutate_conf2(output_type& candidate, output_type& candidate_1, const model& m, fl amplitude, rng& generator, pso* particle, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par,int step, int num_steps, int opt, int count) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion

	output_type tmp_1 = candidate;
	output_type tmp_2 = candidate;

	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	srand((unsigned)time(NULL));
	double a_cm, w_cm, r_cm, alpha_cm, beta_cm, isupdate;
	int y, mc_par=0;

	switch(opt)
	{
		case 1:
			r_cm = (double) rand() / (RAND_MAX + 1.0);
			w_cm = (double) rand() / (RAND_MAX + 1.0);
			r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);
			w_cm = 1.07*(7.86*w_cm-23.31*w_cm*w_cm+28.75*w_cm*w_cm*w_cm-13.302875*w_cm*w_cm*w_cm*w_cm);
			break;
		case 2:
			a_cm = (0.9-0.2)*(1-double(step)/num_steps)+0.2;
			break;
		case 3:
			alpha_cm  = (1.0 - 0.3) * (1-step/(num_steps*0.9)) + 0.3;
			if(alpha_cm < 0.3)
				alpha_cm=0.3;
			beta_cm  = 1.5 -  0.3 * (step/(num_steps*0.9));
			if(beta_cm < 1.2)
				beta_cm = 1.2;
			break;
		default:
			std::cerr << '\n' << "ERROR: opt parameter choose from range of 1 to 6" << '\n';
			VINA_CHECK(false);
	}

	VINA_FOR_IN(i, candidate.c.ligands) {

		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);

		mc_par = random_int(0, particle->number-1, generator);

		for (y=0;y<particle->number;y++)
		{

			candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
			candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);

			isupdate  = (double) rand() / (RAND_MAX + 1.0);
			if(isupdate<0.06)
			{
				quasi_newton_par(tmp_m, p, ig, candidate, g, v);
				particle->updateCurrentPosition(y,candidate.c.ligands[i].rigid.position);
				particle->updateCurrentOrientation(y,candidate.c.ligands[i].rigid.orientation);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					particle->updateCurrentTorsion(y, candidate.c.ligands[i].torsions[z],z);
			}
			else
				candidate.e=tmp_m.eval_deriv(p, ig, v, candidate.c, g);

			particle->updateCurrentBest(y,candidate.e);
			tmp_2 = candidate;
			
			isupdate = (double) rand() / (RAND_MAX + 1.0);
			if(y==mc_par && isupdate<0.01)
			{
				candidate.c.ligands[i].rigid.position =  particle->getPBestPosition(y);
				candidate.c.ligands[i].rigid.orientation = particle->getPBestOrientation(y);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->getPBestTorsion(y,z);
				markov_mutate_conf(candidate, tmp_m, 2, generator, p, ig, g, v, quasi_newton_par);
				if(candidate.e<tmp_2.e)
					tmp_2 = candidate;
			}

			if (tmp_2.e < particle->getPersonalBest(y))
			{
				particle->updatePersonalBest(y,tmp_2.e);
				particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
				particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
				for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
					particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
			}

			isupdate = (double) rand() / (RAND_MAX + 1.0);
			if(y==mc_par && isupdate<0.01)
			{
				candidate.c.ligands[i].rigid.position =  particle->gbest_position;
				candidate.c.ligands[i].rigid.orientation = particle->gbest_orientation;
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->gbest_torsion[z];
				markov_mutate_conf(candidate, tmp_m, 2, generator, p, ig, g, v, quasi_newton_par);
				if(candidate.e<tmp_2.e)
					tmp_2 = candidate;
			}

			if(tmp_2.e < particle->gbest_fit)
			{
				particle->updateGlobalBest_1(tmp_2.e);
				particle->gbest_position = tmp_2.c.ligands[i].rigid.position;
				particle->gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
				for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
					particle->gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
			}

			sz which1 = sz(which_int);
			VINA_CHECK(which1 < mutable_entities_num);

			if(which1 == 0)
			{
				switch(opt)
				{
					case 1:
						particle->updateVelocity(generator,y,r_cm,w_cm);
						particle->computePSOPositions(y);
						break;
					case 2:
						particle->computeQPSOPositions(y,a_cm);
						break;
					case 3:
						particle->computeRDPSOPositions(y,alpha_cm,beta_cm);
						break;
					default:
						std::cerr << '\n' << "ERROR: opt parameter choose from range of 1 to 6" << '\n';
						VINA_CHECK(false);
				}
			}
			--which1;

			if(which1 == 0)
			{
				fl gr = m.gyration_radius(i); 
				if(gr > epsilon_fl)
				{
					switch(opt)
					{
						case 1:
							particle->updateVelocityO(generator,y,r_cm,w_cm);
							particle->computePSOOrientation(y);
							break;
						case 2:
							particle->computeQPSOOrientation(y,a_cm);
							break;
						case 3:
							particle->computeRDPSOOrientation(y,alpha_cm,beta_cm);
							break;
						default:
							std::cerr << '\n' << "ERROR: opt parameter choose from range of 1 to 6" << '\n';
							VINA_CHECK(false);
					}
				}
			}
			--which1;

			if(which1 < candidate.c.ligands[i].torsions.size() && which1>=0)
			{
				switch(opt)
				{
					case 1:
						particle->updateVelocityT(generator,y,which1,r_cm,w_cm);
						particle->computePSOTorsion(y,generator,which1);
						break;
					case 2:
						particle->computeQPSOTorsion(y,which1,a_cm);
						break;
					case 3:
						particle->computeRDPSOTorsion(y,which1,alpha_cm,beta_cm);
						break;
					default:
						std::cerr << '\n' << "ERROR: opt parameter choose from range of 1 to 6" << '\n';
						VINA_CHECK(false);
				}
			}
			which1 -= candidate.c.ligands[i].torsions.size();

		}
	for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
		candidate_1.c.ligands[i].torsions[z] = particle->gbest_torsion[z];
	candidate_1.c.ligands[i].rigid.orientation = particle->gbest_orientation;
	candidate_1.c.ligands[i].rigid.position = particle->gbest_position;
	candidate_1.e = particle->gbest_fit;
	return; 
	}

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}


void pso_mutate_conf1(output_type& candidate, output_type& candidate_1, const model& m, fl amplitude, rng& generator, pso* particle, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par,int step, int num_steps, int count) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion

	output_type tmp_1 = candidate;
	output_type tmp_2 = candidate;

	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	srand((unsigned)time(NULL));
	double  w_cm, r_cm, isupdate;
	int y, mc_par=0;

	r_cm = (double) rand() / (RAND_MAX + 1.0);
	r_cm = 1.07*(7.86*r_cm-23.31*r_cm*r_cm+28.75*r_cm*r_cm*r_cm-13.302875*r_cm*r_cm*r_cm*r_cm);

	w_cm = 0.9 - (double)step / num_steps * 0.5;

	VINA_FOR_IN(i, candidate.c.ligands) {

		model tmp_m = m;
		const vec authentic_v(1000, 1000, 1000);

		mc_par = random_int(0, particle->number-1, generator);

		for (y=0;y<particle->number;y++)
		{

			candidate.c.ligands[i].rigid.position = particle->getCurrentPosition(y);
			candidate.c.ligands[i].rigid.orientation = particle->getCurrentOrientation(y);
			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
				candidate.c.ligands[i].torsions[z] = particle->getCurrentTorsion(y,z);

			isupdate  = (double) rand() / (RAND_MAX + 1.0);
			if(isupdate<0.06)
			{
				quasi_newton_par(tmp_m, p, ig, candidate, g, v);
				particle->updateCurrentPosition(y,candidate.c.ligands[i].rigid.position);
				particle->updateCurrentOrientation(y,candidate.c.ligands[i].rigid.orientation);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					particle->updateCurrentTorsion(y, candidate.c.ligands[i].torsions[z],z);
			}
			else
			{
				if(particle->getPersonalBest(y)>50)
				{
					quasi_newton_par(tmp_m, p, ig, candidate, g, v);
					particle->updateCurrentPosition(y,candidate.c.ligands[i].rigid.position);
					particle->updateCurrentOrientation(y,candidate.c.ligands[i].rigid.orientation);
					for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
						particle->updateCurrentTorsion(y, candidate.c.ligands[i].torsions[z],z);
				}
				else
				{
					candidate.e=tmp_m.eval_deriv(p, ig, v, candidate.c, g);
				}
			}

			particle->updateCurrentBest(y,candidate.e);
			tmp_2 = candidate;

			if(y==mc_par && count<particle->number/2)
			{
				candidate.c.ligands[i].rigid.position =  particle->getPBestPosition(y);
				candidate.c.ligands[i].rigid.orientation = particle->getPBestOrientation(y);
				for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
					candidate.c.ligands[i].torsions[z] = particle->getPBestTorsion(y,z);
				markov_mutate_conf(candidate, tmp_m, 2, generator, p, ig, g, v, quasi_newton_par);
				if(candidate.e<tmp_2.e)
					tmp_2 = candidate;
			}

			if (tmp_2.e < particle->getPersonalBest(y))
			{
				particle->updatePersonalBest(y,tmp_2.e);
				particle->updateBestPosition(y,tmp_2.c.ligands[i].rigid.position);
				particle->updateBestOrientation(y,tmp_2.c.ligands[i].rigid.orientation);
				for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
					particle->updateBestTorsion(y, tmp_2.c.ligands[i].torsions[z],z);
			}

			if(tmp_2.e < particle->gbest_fit)
			{
				particle->updateGlobalBest_1(tmp_2.e);
				particle->gbest_position = tmp_2.c.ligands[i].rigid.position;
				particle->gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
				for(int z=0;z<tmp_2.c.ligands[i].torsions.size();z++)
					particle->gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
			}

			particle->updateVelocity(generator,y,r_cm,w_cm);
			particle->computePSOPositions(y);
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl)
			{
				particle->updateVelocityO(generator,y,r_cm,w_cm);
				particle->computePSOOrientation(y);
			}
			for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
			{
				particle->updateVelocityT(generator,y,z,r_cm,w_cm);
				particle->computePSOTorsion(y,generator,z);
			}

			}

		for(int z=0;z<candidate.c.ligands[i].torsions.size();z++)
			candidate_1.c.ligands[i].torsions[z] = particle->gbest_torsion[z];
		candidate_1.c.ligands[i].rigid.orientation = particle->gbest_orientation;
		candidate_1.c.ligands[i].rigid.position = particle->gbest_position;
		candidate_1.e = particle->gbest_fit;
		return; 
		}

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}

void markov_mutate_conf(output_type& candidate, const model& m, fl amplitude, rng& generator, const precalculate& p ,const igrid& ig,change& g,const vec& v,quasi_newton& quasi_newton_par)
{ // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	output_type tmp_1 = candidate;
	output_type tmp_2 = candidate;
	sz mutable_entities_num = pso_count_mutable_entities(candidate.c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	int which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	VINA_FOR_IN(i, candidate.c.ligands)
	{
		model tmp_m = m;
		if(which == 0){
			candidate.c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator);
		}
		--which;

		if(which == 0)
		{
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl) {
				vec rotation; 
				rotation = amplitude / gr * random_inside_sphere(generator); 
				quaternion_increment(candidate.c.ligands[i].rigid.orientation, rotation);
			}
		}
		--which;

		if(which < candidate.c.ligands[i].torsions.size() && which>=0) 
			candidate.c.ligands[i].torsions[which] = random_fl(-pi, pi, generator);
		
		which -= candidate.c.ligands[i].torsions.size();

		candidate.e=tmp_m.eval_deriv(p, ig, v, candidate.c, g);

		tmp_1 = candidate;
		quasi_newton_par(tmp_m, p, ig, tmp_1, g, v, 1);
		if(tmp_1.e < candidate.e)
		{
			candidate = tmp_1;
			tmp_2=candidate;
			quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);
			if(tmp_2.e < candidate.e)
				candidate = tmp_2;
		}

	}

	VINA_FOR_IN(i, candidate.c.flex) {
		if(which < candidate.c.flex[i].torsions.size()) { candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= candidate.c.flex[i].torsions.size();
	}
}
