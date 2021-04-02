/*
        PSOVina version 1.0                     Date: 26/11/2014

        This file is revised from monte_carlo.cpp in AutoDock Vina.

        Authors: Marcus C. K. Ng  <marcus.ckng@gmail.com>

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


#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "ssd.h"
#include "incrementable.h"
#include "mpso.h"

struct parallel_pso_task {
  model m;
  output_container out;
    int num;
    fl best_fit;
    vec best_position;
    qt best_orientation;
    fl best_torsion[100];
    fl best_fit1;
    vec best_position1;
    qt best_orientation1;
    fl best_torsion1[100];
  rng generator;
  parallel_pso_task(const int num_, const model& m_, int seed) : num(num_), m(m_), generator(static_cast<rng::result_type>(seed)) , best_fit(1.7976931348623158e+308), best_fit1(1.7976931348623158e+308){}
};

typedef boost::ptr_vector<parallel_pso_task> parallel_pso_task_container;
typedef boost::ptr_vector<parallel_pso_task>& parallel_pso_task_container1;

struct pso_search {
	unsigned num_steps;
	fl temperature;
	vec hunt_cap;
	fl min_rmsd;
	sz num_saved_mins;
	fl mutation_amplitude;
	ssd ssd_par;

	pso_search() : num_steps(2500), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2) {} // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  num_steps = 50*lig_atoms = 2500

	output_type operator()(parallel_pso_task_container1 t_vec, int, model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, int num_birds, int opt, double c1, double c2) const;
	// out is sorted
	void operator()(parallel_pso_task_container1 t_vec, int, model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, int num_birds, int opt, double c1, double c2) const;


};

#endif

