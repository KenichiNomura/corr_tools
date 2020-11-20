#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <set>
#include <regex>

#include <cstdlib>
#include <cmath>
#include <cassert>

#include <boost/histogram.hpp>
#include <boost/format.hpp>
using namespace boost::histogram;

namespace fs = std::filesystem;

struct bond_stat
{
	double d0d0 = 0.0;
	std::vector<double> d0dt; 

	std::vector<double> p1, p2; 
};

struct dist 
{
	double dx,dy,dz;
	dist(double _dx, double _dy, double _dz) : dx(_dx), dy(_dy), dz(_dz){};
};

typedef std::tuple<double, int, std::string, dist, bool, int, std::string> tuple_t;
typedef std::vector<std::vector<tuple_t>> nbrlist_t;

#define MAX_NEIGHBOR_MOLS 5
#define HBOND_OUTER_CUTOFF 2.1
#define HBOND_INNER_CUTOFF 1.2

bool nbrlist_compare(nbrlist_t &t1, nbrlist_t &t2)
{
	// return true if two nbrlists are empty 
	if(t1.size()==0 && t2.size()==0) return true;

	// check if number of atoms are the same
	if (t1.size() != t2.size()) return false;

	std::set<int> set1, set2;

	for(int i = 0; i<t1.size(); i++)
	{
		// check if number of nbrs of i-th atom are the same
		if(t1[i].size() != t2[i].size()) return false; 

		// check if the nbr global id is the same
		for(int j = 0; j < t1[i].size(); j++)
		{
			set1.insert(std::get<1>(t1[i][j]));
			set2.insert(std::get<1>(t2[i][j]));
		}
		if(set1 != set2) return false; 
	}

	return true; 
}

// to read .bnd file
typedef std::vector<std::vector<int>> simple_nbrlist_t;

simple_nbrlist_t get_simple_nbrlists(std::string const & filename)
{
	std::ifstream fin(filename);
	simple_nbrlist_t simple_nbrlists; 

	std::string buffer;

	int num_atoms = 0; 
	while(std::getline(fin, buffer)) {num_atoms++;};

	simple_nbrlists.resize(num_atoms);

	fin.clear();
	fin.seekg(0);

	std::stringstream ss;
	while(std::getline(fin, buffer))
	{
		ss << buffer;

		int id,iid,jid,type_,num_nbrs; 
		double x_,y_,z_,dr;
		ss >> id >> x_ >> y_ >> z_ >> type_ >> num_nbrs;

		iid = id-1;
		for(int i=0; i<num_nbrs; i++)
		{
			ss >> jid >> dr; 
			int jid1 = jid - 1; // zero-indexed
			simple_nbrlists[iid].push_back(jid1);
		}

		//std::cout << buffer << std::endl;; 
		ss.str(""); ss.clear();
	}

	return simple_nbrlists;
}

struct MDFrame
{
    std::string filename;

    int natoms;
    double lattice[6];

    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<int> mol_id;  // molecular Id
    std::vector<std::string> name;

    double apply_pbc(double x, double lattice)
    {
        if(x>=0.5*lattice) x-=lattice;
        if(x<-0.5*lattice) x+=lattice;
        return x;
    };

    void print()
    {
      std::cout << " # of Atoms : " << natoms << std::endl;
      std::cout << " Lattice Consts : " <<
              lattice[0] << " " << lattice[1] << " " << lattice[2] << " " <<
              lattice[3] << " " << lattice[4] << " " << lattice[5] << std::endl;
    }
};

struct NeighborList
{
	int max_neighbors;

	nbrlist_t nbrlist;
	nbrlist_t nbrlist_ocenter; 
	nbrlist_t nbrlist_hcenter; 

	simple_nbrlist_t slist;

	NeighborList() {}; 

	NeighborList(MDFrame & mdframe, int _max=MAX_NEIGHBOR_MOLS) : max_neighbors(_max)
	{
		std::string bndfile = std::regex_replace(mdframe.filename, std::regex("xyz"), "bnd");

		if(fs::exists(bndfile)) 
		{
			std::cout << mdframe.filename << " " << bndfile << " " << fs::exists(bndfile) << std::endl;
			slist = get_simple_nbrlists(bndfile); 
		} else {
			std::cout << mdframe.filename << " " << fs::exists(bndfile) << std::endl;
		}

		nbrlist.resize(mdframe.natoms);	
		nbrlist_ocenter.resize(mdframe.natoms);	
		nbrlist_hcenter.resize(mdframe.natoms);	

		auto all_lists = {nbrlist, nbrlist_ocenter, nbrlist_hcenter};

		for (int i=0; i<mdframe.natoms; i++)
		{
			//std::cout << mdframe.name[i] << " " << mdframe.x[i] << " " << mdframe.y[i] << " " << mdframe.z[i] << std::endl; 

			//// construct nbrlist from only O
			//if (mdframe.name[i].find("O") == std::string::npos) continue; 
			//
			bool is_i_oxygen = mdframe.name[i].find("O") != std::string::npos;

			std::vector<int> idlist; 
			if(slist.size()==mdframe.natoms)
			{
				idlist = slist[i];
			} else {
				for(int i=0; i<mdframe.natoms; i++) idlist.push_back(i);
			}

/*
			for (int j1=0; j1<idlist.size(); j1++)
				std::cout << idlist[j1] << " "; 
				std::cout << std::endl;
*/
			for (int j1=0; j1<idlist.size(); j1++)
			{
				int j = idlist[j1];
				if(i==j) continue;

				//std::cout << j << " ";

				bool is_j_oxygen = mdframe.name[j].find("O") != std::string::npos;

				// do not include WCs
				if (mdframe.name[j].find("X") != std::string::npos) continue;

				double dx = mdframe.x[j] - mdframe.x[i];
				double dy = mdframe.y[j] - mdframe.y[i];
				double dz = mdframe.z[j] - mdframe.z[i];

				dx = mdframe.apply_pbc(dx, mdframe.lattice[0]);
				dy = mdframe.apply_pbc(dy, mdframe.lattice[1]);
				dz = mdframe.apply_pbc(dz, mdframe.lattice[2]);

				double dr = sqrt(dx*dx + dy*dy + dz*dz); 
				if (dr == 0.0) continue;

				auto data = std::make_tuple(dr, j, mdframe.name[j], dist(dx,dy,dz), true, i, mdframe.name[i]);

				nbrlist[i].push_back(data);

				// ignore if an atomic pair distance is beyond than h-bond outer cutoff
				if (HBOND_OUTER_CUTOFF < dr) continue;

				if(is_i_oxygen && not is_j_oxygen) nbrlist_ocenter[i].push_back(data);
				if(not is_i_oxygen && is_j_oxygen) nbrlist_hcenter[i].push_back(data);
			}
		}

		// sort all lists by Id
		for (auto list : all_lists)
			for (std::vector<tuple_t> & n : list)
				std::sort(begin(n), end(n), [](auto const & t1, auto const & t2) { return std::get<1>(t1) < std::get<1>(t2);} );

	}

	void print()
	{
		for (int i=0; i<nbrlist.size(); i++)
		{
			std::cout << std::endl << i << "th neighbr list size : " << nbrlist[i].size() << std::endl << std::endl;
			for (const auto & l : nbrlist[i])
				std::cout << i << " " << std::get<0>(l) << " " << std::get<1>(l) << " " << std::get<2>(l) << std::endl; 
				std::cout << std::endl; 
		}
	}

};


MDFrame read_single_mdframe(std::ifstream &in, std::string _filename="NA")
{
    MDFrame mdframe;
    std::string str;

    std::getline(in,str);
    mdframe.filename = _filename;
    mdframe.natoms = std::atoi(str.c_str());

    double dummy;
    std::stringstream ss;
    std::getline(in,str);
    ss << str;
/*
    ss >> str >> mdframe.lattice[0] >> dummy >> dummy >>
      dummy >> mdframe.lattice[1] >> dummy >>
      dummy >> dummy >> mdframe.lattice[2];
*/
    ss >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2];
    mdframe.lattice[3] = mdframe.lattice[4] = mdframe.lattice[5] = 90.0;

    mdframe.name.resize(mdframe.natoms);
    mdframe.x.resize(mdframe.natoms);
    mdframe.y.resize(mdframe.natoms);
    mdframe.z.resize(mdframe.natoms);
    mdframe.mol_id.resize(mdframe.natoms);
    //std::cout << mdframe.natoms << std::endl;

    for (int i=0; i<mdframe.natoms; i++)
    {
        std::string name;
        float x,y,z,vx,vy,vz;
        int id;

        std::stringstream ss;
        std::getline(in,str);
        ss << str;

//        ss >> name >> x >> y >> z >> vx >> vy >> vz;
//
        ss >> name >> x >> y >> z >> dummy >> id;
//        ss >> name >> x >> y >> z >> id; 
        mdframe.name[id-1] = name;
        mdframe.x[id-1] = x;
        mdframe.y[id-1] = y;
        mdframe.z[id-1] = z;
        mdframe.mol_id[id-1] = id;

/*
        ss >> name >> x >> y >> z >> id;
        
        mdframe.name[i] = name;
        mdframe.x[i] = x;
        mdframe.y[i] = y;
        mdframe.z[i] = z;
        mdframe.mol_id[i] = id;
*/

        //std::cout << mdframe.name[i] << " " << mdframe.x[i] << " " << mdframe.y[i] << " " << mdframe.z[i] << std::endl; 
        //std::cout << name << " " <<x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << std::endl;
    }

    return mdframe;
};

auto angle_hoh_ll = make_histogram(axis::regular<double>(180, 0.0, 180.0));
auto angle_hoh_sl = make_histogram(axis::regular<double>(180, 0.0, 180.0));
auto angle_hoh_ss = make_histogram(axis::regular<double>(180, 0.0, 180.0));
auto angle_oho_ll = make_histogram(axis::regular<double>(180, 0.0, 180.0));
auto angle_oho_sl = make_histogram(axis::regular<double>(180, 0.0, 180.0));

class hbond_stats
{
	std::string datadir;
	std::vector<fs::path> filepaths; 

	NeighborList first_nbr; 

	bond_stat bstat; 

	std::string file_histogram, file_correlation;

        // num_alive_bonds, ratio, d0dt, p1, p2
	typedef std::tuple<int, double, double, double, double> corr_t;
	std::map<int,corr_t> corr;

public:

	hbond_stats(std::string _datadir, std::string file_corr, std::string file_histo) 
	: datadir(_datadir), file_correlation(file_corr), file_histogram(file_histo)
	{
		auto cpath = fs::current_path();
		std::cout << cpath.string() << "\n";

		for (auto f : fs::recursive_directory_iterator(datadir))
		{
			if(f.path().string().find(".xyz") != std::string::npos)
			{
				filepaths.push_back(f.path());
				//std::cout << f.path().string() << "\n"; 
			}
		}

		std::sort(begin(filepaths),end(filepaths), 
			[](auto const &p1, auto const &p2) {return p1.string() < p2.string();} );


		// get initial frame
		std::string filename = filepaths[0].string();
		std::ifstream ifile(filename);	
		MDFrame single_frame = read_single_mdframe(ifile, filename);
		first_nbr = NeighborList(single_frame);
	}

	void calc_corr(int interval)
	{
		std::ofstream fout_correlation(file_correlation);	

		for (int step = 0; step<filepaths.size(); step++)
		{
			if(step%interval != 0) continue;
	
			//std::string filename = f.string();
			std::string filename = filepaths[step].string();
			std::ifstream ifile(filename);	
			MDFrame single_frame = read_single_mdframe(ifile, filename);
			//single_frame.print();
	
			auto nbr = NeighborList(single_frame);
	
			assert(first_nbr.nbrlist_ocenter.size() == nbr.nbrlist_ocenter.size()); 
	
			int num_broken_bonds = 0; 
			int num_alive_bonds = 0; 
	
			// bond-bond correlation (single step).  
			double d0dt_step = 0.0;
			int num_d0dt_step = 0;
			double p1_step = 0.0;
			int num_p1_step = 0;
			double p2_step = 0.0;
			int num_p2_step = 0;
	
			std::cout << filename <<  " ===========================\n";

	//=== distribution calculations ===//
			for (int iatom = 0; iatom< nbr.nbrlist_hcenter.size(); iatom++)
			{
				const auto n = nbr.nbrlist_hcenter[iatom]; 
				if(n.size() >= 2) 
				{
					for (int j = 0; j < n.size()-1; j++)
					{
						bool is_j_long = std::get<0>(n[j]) > HBOND_INNER_CUTOFF;

						const auto nj_r = std::get<3>(n[j]);
						const double nj_norm = std::sqrt(nj_r.dx*nj_r.dx + nj_r.dy*nj_r.dy + nj_r.dz*nj_r.dz);
	
						const int jgid = std::get<1>(n[j]);
	
						for(int k = j+1; k < n.size(); k++)
						{
							bool is_k_long = std::get<0>(n[k]) > HBOND_INNER_CUTOFF;

							const auto nk_r = std::get<3>(n[k]);
							const double nk_norm = std::sqrt(nk_r.dx*nk_r.dx + nk_r.dy*nk_r.dy + nk_r.dz*nk_r.dz);
	
							const int kgid = std::get<1>(n[k]);
	
							double dot = nj_r.dx*nk_r.dx + nj_r.dy*nk_r.dy + nj_r.dz*nk_r.dz;
							double cosine = dot/nj_norm/nk_norm;
							double degree = std::acos(cosine) * 180.0 / M_PI;
	
							if(is_j_long && is_k_long) angle_oho_ll(degree);
							if((is_j_long && not is_k_long) || (not is_j_long && is_k_long)) angle_oho_sl(degree);
						}
					}
				}
			}

			for (int iatom = 0; iatom< nbr.nbrlist_ocenter.size(); iatom++)
			{
				const auto n = nbr.nbrlist_ocenter[iatom]; 
				if(n.size() >= 2) 
				{
					for (int j = 0; j < n.size()-1; j++)
					{
						bool is_j_long = std::get<0>(n[j]) > HBOND_INNER_CUTOFF;

						const auto nj_r = std::get<3>(n[j]);
						const double nj_norm = std::sqrt(nj_r.dx*nj_r.dx + nj_r.dy*nj_r.dy + nj_r.dz*nj_r.dz);
	
						const int jgid = std::get<1>(n[j]);

	
						for(int k = j+1; k < n.size(); k++)
						{
							bool is_k_long = std::get<0>(n[k]) > HBOND_INNER_CUTOFF;

							const auto nk_r = std::get<3>(n[k]);
							const double nk_norm = std::sqrt(nk_r.dx*nk_r.dx + nk_r.dy*nk_r.dy + nk_r.dz*nk_r.dz);
	
							const int kgid = std::get<1>(n[k]);
	
							double dot = nj_r.dx*nk_r.dx + nj_r.dy*nk_r.dy + nj_r.dz*nk_r.dz;
							double cosine = dot/nj_norm/nk_norm;
							double degree = std::acos(cosine) * 180.0 / M_PI;
	
							if(is_j_long && is_k_long) angle_hoh_ll(degree);
							if(not is_j_long && not is_k_long) angle_hoh_ss(degree);
							if((is_j_long && not is_k_long) || (not is_j_long && is_k_long)) angle_hoh_sl(degree);
						}
					}
				}
			}

	
	//=== correlation calculations ===//
			//for (auto n : nbr.nbrlist_ocenter)
			for (int iatom = 0; iatom< nbr.nbrlist_ocenter.size(); iatom++)
			{
				auto m = first_nbr.nbrlist_ocenter[iatom]; 
				const auto n = nbr.nbrlist_ocenter[iatom]; 
	
	
				if(m.size() > 0 && n.size() > 0) 
				{
					//std::cout << m.size() << " " << n.size() << " " << "====================\n";
	
					for (int j = 0; j < m.size(); j++)
					{
						// this bond is already broken, do not sample it.
						if(not std::get<4>(m[j])) continue; 
	
						const auto m_r = std::get<3>(m[j]);
						const double m_norm = std::sqrt(m_r.dx*m_r.dx + m_r.dy*m_r.dy + m_r.dz*m_r.dz);
	
						const int jgid = std::get<1>(m[j]);
	
						// first set the bond alive flag false, then reset to true if it is bonded.
						std::get<4>(m[j]) = false; 
	
						for(int k = 0; k < n.size(); k++)
						{
							const auto n_r = std::get<3>(n[k]);
							const double n_norm = std::sqrt(n_r.dx*n_r.dx + n_r.dy*n_r.dy + n_r.dz*n_r.dz);
	
							const int kgid = std::get<1>(n[k]);
	
							if(jgid == kgid) // compute bond-correlation
							{
								double dot = m_r.dx*n_r.dx + m_r.dy*n_r.dy + m_r.dz*n_r.dz;
								double cosine = dot/m_norm/n_norm;
	
	/*
								std::cout << std::get<1>(m[j]) << " " << std::get<1>(n[k]) << " " << 
									m_r.dx << " " << n_r.dx << " " << 
									m_r.dy << " " << n_r.dy << " " <<  
									m_r.dz << " " << n_r.dz << " " << dot << std::endl;
	*/
	
								d0dt_step += dot; 
								num_d0dt_step++; 
	
								std::get<4>(m[j]) = true; 
	
								p1_step += cosine; 
								num_p1_step++;
								p2_step += 0.5*(3*std::pow(cosine,2)-1.0);
								num_p2_step++;
							}
						}
	
						if(std::get<4>(m[j]))
							num_alive_bonds++;  
						else 
							num_broken_bonds++;  
					}
				}
			}
	
			bstat.d0dt.push_back(d0dt_step/num_d0dt_step); 
			bstat.p1.push_back(p1_step/num_p1_step); 
			bstat.p2.push_back(p2_step/num_p2_step);
	
			if(step==0) bstat.d0d0 = bstat.d0dt.back(); 
	
			fout_correlation << "step,num_alive_bonds,ratio,d0dt,p1,p2 : " << step << " " << 
				num_alive_bonds << " " << (double)num_alive_bonds/(num_alive_bonds+num_broken_bonds) << " " << 
				bstat.d0dt.back()/bstat.d0d0 << " " << bstat.p1.back() << " " << bstat.p2.back() << std::endl;

				
			corr[step] = std::make_tuple(num_alive_bonds, 
				(double)num_alive_bonds/(num_alive_bonds+num_broken_bonds), 
				bstat.d0dt.back()/bstat.d0d0, bstat.p1.back(), bstat.p2.back() );
	
		}

		fout_correlation.close();
	};

	void save_result()
	{
		std::vector<std::vector<double>> angles;

		// normalize all histograms
		auto all_histograms = {angle_hoh_ll, angle_hoh_sl, angle_hoh_ss, angle_oho_ll, angle_oho_sl};
		for(auto h : all_histograms)
		{
			std::vector<double> angle;

			double sum = std::accumulate(h.begin(), h.end(), 0.0);
			for(int i=0; i<h.axis().size(); i++) angle.push_back(h[i]/sum);
			angles.push_back(angle);
                }

		std::ofstream fout_histogram(file_histogram);	
	
		fout_histogram << "angle    HOH_ll    HOH_sl    HOH_ss    OHO_ll    OHO_sl" << std::endl;
		for(int idx=0; idx<angles[0].size(); idx++) 
		{
			fout_histogram << idx << " ";
			for(int j=0; j<angles.size(); j++) fout_histogram << angles[j][idx] << " ";
			fout_histogram << std::endl;
		}

		fout_histogram.close();

		std::ofstream fout_correlation(file_correlation);	
		fout_correlation << "step   bond_avlie  bond_ratio  d0dt  p1  p2 " << std::endl;
		for(auto &[step,c] : corr)
			fout_correlation << step << " " << std::get<0>(c) << " " << std::get<1>(c) << " " << 
				std::get<2>(c) << " " << std::get<3>(c) << " " << std::get<4>(c) << " " << std::endl;

		fout_correlation.close();
	};

};

int main(int argc, char* argv[])
{
	std::string datadir(argv[1]);
	int interval = std::atoi(argv[2]);

	std::cout << "datadir, interval : " << datadir << " " << interval << std::endl << std::endl;

	hbond_stats hstats(datadir, "corr.dat", "hbond.dat");
	hstats.calc_corr(interval);
	hstats.save_result();

	return 0;
}
